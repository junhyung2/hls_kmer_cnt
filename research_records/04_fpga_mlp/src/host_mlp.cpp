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

struct ht_entry_mlp_t {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
};

#ifndef HT_BITS_HOST_MLP
  #define HT_BITS_HOST_MLP 22
#endif

constexpr uint64_t HT_T_SIZE  = 1ULL << (HT_BITS_HOST_MLP - 1);   /* T_SIZE per CU */
constexpr uint64_t HT_T_BYTES = HT_T_SIZE * sizeof(ht_entry_mlp_t);
static const int   N_CU       = 4;

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

    for (int i = 5; i < argc - 1; i++)
        if (std::string(argv[i]) == "-d") device_idx = std::stoul(argv[i+1]);

    if (k < 1 || k > 31) { std::cerr << "k must be in [1,31]\n"; return 1; }

    std::cout << "[kmer_count_mlp] Loading FASTQ files ...\n";
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    uint64_t len1 = seq1.size(), len2 = seq2.size();
    std::cout << "  FASTQ1 : " << fq1_path << "  (" << len1 << " B)\n";
    std::cout << "  FASTQ2 : " << fq2_path << "  (" << len2 << " B)\n";
    std::cout << "  k      : " << k << "\n";
    std::cout << "  N_CU   : " << N_CU << "\n";
    std::cout << "  HT_BITS: " << HT_BITS_HOST_MLP
              << "  (T_SIZE=" << HT_T_SIZE << " entries = "
              << HT_T_BYTES/(1024*1024) << " MB per T per CU)\n";
    std::cout << "  Total HT: " << (N_CU * 2 * HT_T_BYTES)/(1024*1024) << " MB\n";

    std::cout << "[kmer_count_mlp] Opening device " << device_idx << " ...\n";
    auto device = xrt::device(device_idx);
    auto uuid   = device.load_xclbin(xclbin_path);

    auto krnl1 = xrt::kernel(device, uuid, "kmer_count_mlp:{kmer_count_mlp_1}");
    auto krnl2 = xrt::kernel(device, uuid, "kmer_count_mlp:{kmer_count_mlp_2}");
    auto krnl3 = xrt::kernel(device, uuid, "kmer_count_mlp:{kmer_count_mlp_3}");
    auto krnl4 = xrt::kernel(device, uuid, "kmer_count_mlp:{kmer_count_mlp_4}");

    /* seq1/seq2: shared read-only across all CUs */
    auto bo_seq1 = xrt::bo(device, len1, krnl1.group_id(0));
    auto bo_seq2 = xrt::bo(device, len2, krnl1.group_id(1));

    /* Per-CU T1 and T2 buffers (kernel args 2 and 3) */
    auto bo_t1_1 = xrt::bo(device, HT_T_BYTES, krnl1.group_id(2));
    auto bo_t2_1 = xrt::bo(device, HT_T_BYTES, krnl1.group_id(3));
    auto bo_t1_2 = xrt::bo(device, HT_T_BYTES, krnl2.group_id(2));
    auto bo_t2_2 = xrt::bo(device, HT_T_BYTES, krnl2.group_id(3));
    auto bo_t1_3 = xrt::bo(device, HT_T_BYTES, krnl3.group_id(2));
    auto bo_t2_3 = xrt::bo(device, HT_T_BYTES, krnl3.group_id(3));
    auto bo_t1_4 = xrt::bo(device, HT_T_BYTES, krnl4.group_id(2));
    auto bo_t2_4 = xrt::bo(device, HT_T_BYTES, krnl4.group_id(3));

    /* Per-CU n_unique outputs (kernel arg 7) */
    auto bo_nu1 = xrt::bo(device, sizeof(uint64_t), krnl1.group_id(7));
    auto bo_nu2 = xrt::bo(device, sizeof(uint64_t), krnl2.group_id(7));
    auto bo_nu3 = xrt::bo(device, sizeof(uint64_t), krnl3.group_id(7));
    auto bo_nu4 = xrt::bo(device, sizeof(uint64_t), krnl4.group_id(7));

    std::memcpy(bo_seq1.map<uint8_t*>(), seq1.data(), len1);
    std::memcpy(bo_seq2.map<uint8_t*>(), seq2.data(), len2);

    std::cout << "[kmer_count_mlp] Initializing hash tables ("
              << (N_CU * 2 * HT_T_BYTES)/(1024*1024) << " MB total) ...\n";

    std::memset(bo_t1_1.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t2_1.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t1_2.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t2_2.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t1_3.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t2_3.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t1_4.map<uint8_t*>(), 0xFF, HT_T_BYTES);
    std::memset(bo_t2_4.map<uint8_t*>(), 0xFF, HT_T_BYTES);

    *bo_nu1.map<uint64_t*>() = 0ULL; *bo_nu2.map<uint64_t*>() = 0ULL;
    *bo_nu3.map<uint64_t*>() = 0ULL; *bo_nu4.map<uint64_t*>() = 0ULL;

    std::cout << "[kmer_count_mlp] Syncing buffers to device ...\n";
    bo_seq1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_seq2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_t1_1.sync(XCL_BO_SYNC_BO_TO_DEVICE); bo_t2_1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_t1_2.sync(XCL_BO_SYNC_BO_TO_DEVICE); bo_t2_2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_t1_3.sync(XCL_BO_SYNC_BO_TO_DEVICE); bo_t2_3.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_t1_4.sync(XCL_BO_SYNC_BO_TO_DEVICE); bo_t2_4.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu1.sync(XCL_BO_SYNC_BO_TO_DEVICE); bo_nu2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu3.sync(XCL_BO_SYNC_BO_TO_DEVICE); bo_nu4.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    /* Launch all 4 CUs simultaneously (args: seq1,seq2,ht_t1,ht_t2,len1,len2,k,n_unique,cu_id) */
    std::cout << "[kmer_count_mlp] Launching 4 CUs ...\n";
    auto t0   = std::chrono::high_resolution_clock::now();
    auto run1 = krnl1(bo_seq1, bo_seq2, bo_t1_1, bo_t2_1, len1, len2, k, bo_nu1, 0);
    auto run2 = krnl2(bo_seq1, bo_seq2, bo_t1_2, bo_t2_2, len1, len2, k, bo_nu2, 1);
    auto run3 = krnl3(bo_seq1, bo_seq2, bo_t1_3, bo_t2_3, len1, len2, k, bo_nu3, 2);
    auto run4 = krnl4(bo_seq1, bo_seq2, bo_t1_4, bo_t2_4, len1, len2, k, bo_nu4, 3);

    run1.wait(); run2.wait(); run3.wait(); run4.wait();
    auto t1 = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[kmer_count_mlp] All CUs finished in " << elapsed << " s\n";

    bo_nu1.sync(XCL_BO_SYNC_BO_FROM_DEVICE); bo_nu2.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu3.sync(XCL_BO_SYNC_BO_FROM_DEVICE); bo_nu4.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    uint64_t nu[4] = {
        *bo_nu1.map<uint64_t*>(), *bo_nu2.map<uint64_t*>(),
        *bo_nu3.map<uint64_t*>(), *bo_nu4.map<uint64_t*>()
    };
    uint64_t total_unique = 0;
    for (int i = 0; i < 4; i++) total_unique += nu[i];

    std::cout << "\n=== Results ===\n";
    for (int i = 0; i < 4; i++)
        std::cout << "  CU" << i << " unique " << k << "-mers : " << nu[i] << "\n";
    std::cout << "  Total unique " << k << "-mers : " << total_unique << "\n";
    std::cout << "  Throughput : "
              << (double)(len1 + len2) / elapsed / 1e6 << " MB/s\n";

    return 0;
}
