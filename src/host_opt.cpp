#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <stdexcept>
#include <cstdlib>

#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

/* Must match kernel header */
struct ht_entry_opt_t {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
};

#ifndef HT_BITS_HOST_OPT
  #define HT_BITS_HOST_OPT 22   /* default: HW_EMU */
#endif
#define N_CU 4

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

    std::cout << "[kmer_count_opt] Loading FASTQ files ...\n";
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    uint64_t len1 = seq1.size();
    uint64_t len2 = seq2.size();
    std::cout << "  FASTQ1 : " << fq1_path << "  (" << len1 << " B)\n";
    std::cout << "  FASTQ2 : " << fq2_path << "  (" << len2 << " B)\n";
    std::cout << "  k      : " << k << "\n";
    std::cout << "  N_CU   : " << N_CU << "\n";

    constexpr uint64_t HT_BITS  = HT_BITS_HOST_OPT;
    constexpr uint64_t HT_SIZE  = 1ULL << HT_BITS;
    constexpr uint64_t HT_BYTES = HT_SIZE * sizeof(ht_entry_opt_t);
    std::cout << "  HT/CU  : 2^" << HT_BITS << " entries = "
              << HT_BYTES / (1024*1024) << " MB\n";

    std::cout << "[kmer_count_opt] Opening device " << device_idx << " ...\n";
    auto device = xrt::device(device_idx);
    auto uuid   = device.load_xclbin(xclbin_path);

    /* Open all 4 CU instances */
    auto krnl1 = xrt::kernel(device, uuid, "kmer_count_opt:{kmer_count_opt_1}");
    auto krnl2 = xrt::kernel(device, uuid, "kmer_count_opt:{kmer_count_opt_2}");
    auto krnl3 = xrt::kernel(device, uuid, "kmer_count_opt:{kmer_count_opt_3}");
    auto krnl4 = xrt::kernel(device, uuid, "kmer_count_opt:{kmer_count_opt_4}");

    /* seq1/seq2 shared read-only buffers (group_id uses krnl1 → HBM[0]/[1]) */
    auto bo_seq1 = xrt::bo(device, len1,    krnl1.group_id(0));
    auto bo_seq2 = xrt::bo(device, len2,    krnl1.group_id(1));

    /* Per-CU hash tables (each routed to its own HBM bank via connectivity.cfg) */
    auto bo_ht1 = xrt::bo(device, HT_BYTES, krnl1.group_id(2));
    auto bo_ht2 = xrt::bo(device, HT_BYTES, krnl2.group_id(2));
    auto bo_ht3 = xrt::bo(device, HT_BYTES, krnl3.group_id(2));
    auto bo_ht4 = xrt::bo(device, HT_BYTES, krnl4.group_id(2));

    /* Per-CU n_unique outputs */
    auto bo_nu1 = xrt::bo(device, sizeof(uint64_t), krnl1.group_id(6));
    auto bo_nu2 = xrt::bo(device, sizeof(uint64_t), krnl2.group_id(6));
    auto bo_nu3 = xrt::bo(device, sizeof(uint64_t), krnl3.group_id(6));
    auto bo_nu4 = xrt::bo(device, sizeof(uint64_t), krnl4.group_id(6));

    /* Copy FASTQ */
    std::memcpy(bo_seq1.map<uint8_t*>(), seq1.data(), len1);
    std::memcpy(bo_seq2.map<uint8_t*>(), seq2.data(), len2);

    /* Initialize hash tables: all EMPTY (0xFF) */
    std::cout << "[kmer_count_opt] Initializing 4x hash tables ("
              << (4 * HT_BYTES) / (1024*1024) << " MB total) ...\n";
    std::memset(bo_ht1.map<uint8_t*>(), 0xFF, HT_BYTES);
    std::memset(bo_ht2.map<uint8_t*>(), 0xFF, HT_BYTES);
    std::memset(bo_ht3.map<uint8_t*>(), 0xFF, HT_BYTES);
    std::memset(bo_ht4.map<uint8_t*>(), 0xFF, HT_BYTES);

    *bo_nu1.map<uint64_t*>() = 0ULL;
    *bo_nu2.map<uint64_t*>() = 0ULL;
    *bo_nu3.map<uint64_t*>() = 0ULL;
    *bo_nu4.map<uint64_t*>() = 0ULL;

    /* Sync all buffers to device */
    std::cout << "[kmer_count_opt] Syncing buffers to device ...\n";
    bo_seq1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_seq2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_ht1 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_ht2 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_ht3 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_ht4 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu1 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu2 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu3 .sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu4 .sync(XCL_BO_SYNC_BO_TO_DEVICE);

    /* Launch all 4 CUs simultaneously */
    std::cout << "[kmer_count_opt] Launching 4 CUs ...\n";
    auto t0   = std::chrono::high_resolution_clock::now();
    auto run1 = krnl1(bo_seq1, bo_seq2, bo_ht1, len1, len2, k, bo_nu1, 0);
    auto run2 = krnl2(bo_seq1, bo_seq2, bo_ht2, len1, len2, k, bo_nu2, 1);
    auto run3 = krnl3(bo_seq1, bo_seq2, bo_ht3, len1, len2, k, bo_nu3, 2);
    auto run4 = krnl4(bo_seq1, bo_seq2, bo_ht4, len1, len2, k, bo_nu4, 3);

    run1.wait();
    run2.wait();
    run3.wait();
    run4.wait();
    auto t1   = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[kmer_count_opt] All CUs finished in " << elapsed << " s\n";

    /* Read back n_unique from each CU */
    bo_nu1.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu2.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu3.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu4.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    uint64_t nu1 = *bo_nu1.map<uint64_t*>();
    uint64_t nu2 = *bo_nu2.map<uint64_t*>();
    uint64_t nu3 = *bo_nu3.map<uint64_t*>();
    uint64_t nu4 = *bo_nu4.map<uint64_t*>();
    uint64_t n_unique_total = nu1 + nu2 + nu3 + nu4;

    std::cout << "\n=== Results ===\n";
    std::cout << "  CU0 unique " << k << "-mers : " << nu1 << "\n";
    std::cout << "  CU1 unique " << k << "-mers : " << nu2 << "\n";
    std::cout << "  CU2 unique " << k << "-mers : " << nu3 << "\n";
    std::cout << "  CU3 unique " << k << "-mers : " << nu4 << "\n";
    std::cout << "  Total unique " << k << "-mers : " << n_unique_total << "\n";
    std::cout << "  Throughput : "
              << (double)(len1 + len2) / elapsed / 1e6 << " MB/s\n";

    /* Optional: dump first 20 occupied entries from CU0 */
    if (std::getenv("KMER_DUMP")) {
        bo_ht1.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        auto* ht = bo_ht1.map<ht_entry_opt_t*>();
        std::cout << "\n[CU0 top-20 occupied entries]\n";
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
