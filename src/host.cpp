#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <stdexcept>
#include <cstdlib>

// XRT native API (Vitis 2022.1 / XRT 2.13)
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

// ht_entry_t definition (mirrors kernel struct)
struct ht_entry_t {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
};

#ifndef HT_BITS_HOST
  #define HT_BITS_HOST 22   // default: SW_EMU
#endif

static std::vector<uint8_t> read_file(const std::string& path)
{
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    auto size = static_cast<size_t>(f.tellg());
    f.seekg(0);
    std::vector<uint8_t> buf(size);
    f.read(reinterpret_cast<char*>(buf.data()), size);
    if (!f) throw std::runtime_error("Read error: " + path);
    return buf;
}

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <xclbin> <fastq1> <fastq2> <k> [-d <device_index>]\n";
        return 1;
    }

    std::string xclbin_path = argv[1];
    std::string fq1_path    = argv[2];
    std::string fq2_path    = argv[3];
    int         k           = std::stoi(argv[4]);
    unsigned    device_idx  = 0;

    for (int i = 5; i < argc - 1; i++) {
        if (std::string(argv[i]) == "-d")
            device_idx = std::stoul(argv[i + 1]);
    }

    if (k < 1 || k > 31) {
        std::cerr << "k must be in [1, 31]\n";
        return 1;
    }

    // ── Load FASTQ ──────────────────────────────────────────────────────
    std::cout << "[kmer_count] Loading FASTQ files ...\n";
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    uint64_t len1 = seq1.size();
    uint64_t len2 = seq2.size();
    std::cout << "  FASTQ1 : " << fq1_path << "  (" << len1 << " B)\n";
    std::cout << "  FASTQ2 : " << fq2_path << "  (" << len2 << " B)\n";
    std::cout << "  k      : " << k << "\n";

    // ── Hash table sizing ────────────────────────────────────────────────
    constexpr uint64_t HT_BITS  = HT_BITS_HOST;
    constexpr uint64_t HT_SIZE  = 1ULL << HT_BITS;
    constexpr uint64_t HT_BYTES = HT_SIZE * sizeof(ht_entry_t);  // 16 B/entry
    std::cout << "  HT     : 2^" << HT_BITS << " entries = "
              << HT_BYTES / (1024*1024) << " MB\n";

    // ── Open FPGA device and load xclbin ─────────────────────────────────
    std::cout << "[kmer_count] Opening device " << device_idx << " and loading xclbin ...\n";
    auto device = xrt::device(device_idx);
    auto uuid   = device.load_xclbin(xclbin_path);
    auto krnl   = xrt::kernel(device, uuid, "kmer_count");

    // ── Allocate XRT buffers (group_id ties to sp= in connectivity.cfg) ──
    // arg 0 → seq1 (gmem0 → HBM[0])
    // arg 1 → seq2 (gmem1 → HBM[1])
    // arg 2 → ht   (gmem2 → HBM[2] or HBM[2:17] for HW)
    // arg 6 → n_unique (gmem3 → HBM[18] or HBM[2] for SW/HW_EMU)
    auto bo_seq1 = xrt::bo(device, len1,    krnl.group_id(0));
    auto bo_seq2 = xrt::bo(device, len2,    krnl.group_id(1));
    auto bo_ht   = xrt::bo(device, HT_BYTES, krnl.group_id(2));
    auto bo_nuni = xrt::bo(device, sizeof(uint64_t), krnl.group_id(6));

    // ── Copy FASTQ data ──────────────────────────────────────────────────
    auto* p1 = bo_seq1.map<uint8_t*>();
    auto* p2 = bo_seq2.map<uint8_t*>();
    std::memcpy(p1, seq1.data(), len1);
    std::memcpy(p2, seq2.data(), len2);

    // ── Initialize hash table: all keys = EMPTY (0xFF...FF) ─────────────
    // memset(0xFF): key=0xFFFF...FFFF (EMPTY_KEY), count/pad=0xFFFFFFFF
    // Count field of empty slots is never read, so this is safe.
    std::cout << "[kmer_count] Initializing hash table ("
              << HT_BYTES / (1024*1024) << " MB) ...\n";
    auto* ht_ptr = bo_ht.map<uint8_t*>();
    std::memset(ht_ptr, 0xFF, HT_BYTES);

    // ── Initialize n_unique output to 0 ─────────────────────────────────
    *bo_nuni.map<uint64_t*>() = 0ULL;

    // ── Sync to device ───────────────────────────────────────────────────
    std::cout << "[kmer_count] Syncing buffers to device ...\n";
    bo_seq1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_seq2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_ht  .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nuni.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // ── Launch kernel ────────────────────────────────────────────────────
    std::cout << "[kmer_count] Launching kernel ...\n";
    auto t0  = std::chrono::high_resolution_clock::now();
    auto run = krnl(bo_seq1, bo_seq2, bo_ht, len1, len2, k, bo_nuni);
    run.wait();
    auto t1  = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[kmer_count] Kernel finished in " << elapsed << " s\n";

    // ── Read back result ─────────────────────────────────────────────────
    bo_nuni.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    uint64_t n_unique = *bo_nuni.map<uint64_t*>();
    std::cout << "\n=== Results ===\n";
    std::cout << "  Unique " << k << "-mers : " << n_unique << "\n";
    std::cout << "  Throughput     : "
              << (double)(len1 + len2) / elapsed / 1e6 << " MB/s\n";

    // ── Optional: dump first 20 occupied entries ─────────────────────────
    bool dump = (std::getenv("KMER_DUMP") != nullptr);
    if (dump) {
        bo_ht.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        auto* ht = bo_ht.map<ht_entry_t*>();
        std::cout << "\n[top-20 occupied entries]\n";
        int shown = 0;
        for (uint64_t i = 0; i < HT_SIZE && shown < 20; i++) {
            if (ht[i].key != 0xFFFFFFFFFFFFFFFFULL) {
                std::cout << "  bucket=" << i
                          << "  key=0x" << std::hex << ht[i].key << std::dec
                          << "  count=" << ht[i].count << "\n";
                shown++;
            }
        }
    }

    return 0;
}
