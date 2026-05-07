#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <stdexcept>
#include <thread>

#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

/* ---- k-mer extraction (CPU) ---- */

static const uint8_t NT4[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask;
    key =   key ^ (key >> 24);
    key = ((key + (key <<  3)) + (key <<  8)) & mask;
    key =   key ^ (key >> 14);
    key = ((key + (key <<  2)) + (key <<  4)) & mask;
    key =   key ^ (key >> 28);
    key =  (key + (key << 31)) & mask;
    return key;
}

/* Extract canonical k-mer hashes from a FASTQ byte buffer.
 * Routes each hash h to bucket (h & 3) matching cu_buckets index. */
static void extract_kmers_fastq(
    const uint8_t*                      seq,
    uint64_t                            len,
    int                                 k,
    uint64_t                            km_mask,
    int                                 shift,
    std::vector<std::vector<uint64_t>>& cu_buckets  /* [4] */
) {
    uint64_t x0 = 0, x1 = 0;
    int l = 0, rline = 0;

    for (uint64_t i = 0; i < len; i++) {
        uint8_t raw = seq[i];
        if (raw == '\n') {
            if (rline == 1) { l = 0; x0 = 0; x1 = 0; }
            rline = (rline == 3) ? 0 : rline + 1;
        } else if (rline == 1) {
            uint8_t c = NT4[raw];
            if (c < 4) {
                x0 = (x0 << 2 | (uint64_t)c) & km_mask;
                x1 = (x1 >> 2) | ((uint64_t)(3 - c) << shift);
                if (++l >= k) {
                    uint64_t y = (x0 < x1) ? x0 : x1;
                    uint64_t h = hash64(y, km_mask);
                    cu_buckets[h & 3].push_back(h);
                }
            } else { l = 0; x0 = 0; x1 = 0; }
        }
    }
}

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

/* ---- HT params (must match kernel) ---- */
#ifndef HT_BITS_HOST_OPT1
  #define HT_BITS_HOST_OPT1 22
#endif
constexpr uint64_t HT_T_SIZE  = 1ULL << (HT_BITS_HOST_OPT1 - 1);
constexpr uint64_t HT_T_BYTES = HT_T_SIZE * 16ULL;  /* sizeof(ht_entry_opt1_t) */
static const int   N_CU       = 4;

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

    uint64_t km_mask = (1ULL << (2*k)) - 1;
    int      shift   = (k-1)*2;

    /* ---- CPU: read FASTQ and extract k-mers ---- */
    std::cout << "[opt1] Loading FASTQ files ...\n";
    auto t_load0 = std::chrono::high_resolution_clock::now();
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    std::cout << "  FASTQ1: " << fq1_path << "  (" << seq1.size() << " B)\n";
    std::cout << "  FASTQ2: " << fq2_path << "  (" << seq2.size() << " B)\n";

    std::cout << "[opt1] CPU k-mer extraction (k=" << k << ") ...\n";
    std::vector<std::vector<uint64_t>> cu_buckets(N_CU);
    extract_kmers_fastq(seq1.data(), seq1.size(), k, km_mask, shift, cu_buckets);
    extract_kmers_fastq(seq2.data(), seq2.size(), k, km_mask, shift, cu_buckets);
    seq1.clear(); seq1.shrink_to_fit();
    seq2.clear(); seq2.shrink_to_fit();

    auto t_load1 = std::chrono::high_resolution_clock::now();
    double t_cpu = std::chrono::duration<double>(t_load1 - t_load0).count();
    uint64_t total_kmers = 0;
    for (int i = 0; i < N_CU; i++) {
        std::cout << "  CU" << i << " k-mers: " << cu_buckets[i].size() << "\n";
        total_kmers += cu_buckets[i].size();
    }
    std::cout << "  Total k-mers: " << total_kmers
              << "  (CPU extraction: " << t_cpu << " s)\n";

    /* ---- XRT setup ---- */
    std::cout << "[opt1] Opening device " << device_idx << " ...\n";
    auto device = xrt::device(device_idx);
    auto uuid   = device.load_xclbin(xclbin_path);

    auto krnl1 = xrt::kernel(device, uuid, "kmer_count_opt1:{kmer_count_opt1_1}");
    auto krnl2 = xrt::kernel(device, uuid, "kmer_count_opt1:{kmer_count_opt1_2}");
    auto krnl3 = xrt::kernel(device, uuid, "kmer_count_opt1:{kmer_count_opt1_3}");
    auto krnl4 = xrt::kernel(device, uuid, "kmer_count_opt1:{kmer_count_opt1_4}");

    /* ---- Allocate k-mer buffers (sized to actual count, no wasted HBM) ---- */
    auto make_kmer_bo = [&](xrt::kernel& krnl, int cu) {
        uint64_t nbytes = cu_buckets[cu].size() * sizeof(uint64_t);
        if (nbytes == 0) nbytes = 8;  /* XRT requires nonzero */
        auto bo = xrt::bo(device, nbytes, krnl.group_id(0));
        if (!cu_buckets[cu].empty())
            std::memcpy(bo.map<uint64_t*>(), cu_buckets[cu].data(), nbytes);
        return bo;
    };

    std::cout << "[opt1] Allocating and syncing k-mer buffers ...\n";
    auto bo_km1 = make_kmer_bo(krnl1, 0);
    auto bo_km2 = make_kmer_bo(krnl2, 1);
    auto bo_km3 = make_kmer_bo(krnl3, 2);
    auto bo_km4 = make_kmer_bo(krnl4, 3);

    /* ---- Allocate HT buffers (NO host init, NO sync to device) ---- */
    std::cout << "[opt1] Allocating HT buffers ("
              << (N_CU * 2 * HT_T_BYTES / 1024 / 1024) << " MB, on-device init) ...\n";
    auto bo_t1_1 = xrt::bo(device, HT_T_BYTES, krnl1.group_id(1));
    auto bo_t2_1 = xrt::bo(device, HT_T_BYTES, krnl1.group_id(2));
    auto bo_t1_2 = xrt::bo(device, HT_T_BYTES, krnl2.group_id(1));
    auto bo_t2_2 = xrt::bo(device, HT_T_BYTES, krnl2.group_id(2));
    auto bo_t1_3 = xrt::bo(device, HT_T_BYTES, krnl3.group_id(1));
    auto bo_t2_3 = xrt::bo(device, HT_T_BYTES, krnl3.group_id(2));
    auto bo_t1_4 = xrt::bo(device, HT_T_BYTES, krnl4.group_id(1));
    auto bo_t2_4 = xrt::bo(device, HT_T_BYTES, krnl4.group_id(2));

    /* n_unique output buffers */
    auto bo_nu1 = xrt::bo(device, sizeof(uint64_t), krnl1.group_id(4));
    auto bo_nu2 = xrt::bo(device, sizeof(uint64_t), krnl2.group_id(4));
    auto bo_nu3 = xrt::bo(device, sizeof(uint64_t), krnl3.group_id(4));
    auto bo_nu4 = xrt::bo(device, sizeof(uint64_t), krnl4.group_id(4));
    *bo_nu1.map<uint64_t*>() = 0; *bo_nu2.map<uint64_t*>() = 0;
    *bo_nu3.map<uint64_t*>() = 0; *bo_nu4.map<uint64_t*>() = 0;

    /* Sync: only k-mer arrays and n_unique init to device (NOT HT) */
    std::cout << "[opt1] Syncing k-mer buffers to device ...\n";
    auto t_sync0 = std::chrono::high_resolution_clock::now();
    bo_km1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_km2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_km3.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_km4.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu3.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nu4.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    auto t_sync1 = std::chrono::high_resolution_clock::now();
    double t_sync = std::chrono::duration<double>(t_sync1 - t_sync0).count();
    std::cout << "  k-mer sync: " << t_sync << " s\n";

    /* ---- Launch all 4 CUs (args: kmers, ht_t1, ht_t2, n_kmers, n_unique) ---- */
    std::cout << "[opt1] Launching 4 CUs (on-device HT init + insert) ...\n";
    auto t0 = std::chrono::high_resolution_clock::now();
    auto run1 = krnl1(bo_km1, bo_t1_1, bo_t2_1, (uint64_t)cu_buckets[0].size(), bo_nu1);
    auto run2 = krnl2(bo_km2, bo_t1_2, bo_t2_2, (uint64_t)cu_buckets[1].size(), bo_nu2);
    auto run3 = krnl3(bo_km3, bo_t1_3, bo_t2_3, (uint64_t)cu_buckets[2].size(), bo_nu3);
    auto run4 = krnl4(bo_km4, bo_t1_4, bo_t2_4, (uint64_t)cu_buckets[3].size(), bo_nu4);
    run1.wait(); run2.wait(); run3.wait(); run4.wait();
    auto t1 = std::chrono::high_resolution_clock::now();
    double t_kernel = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[opt1] Kernel finished in " << t_kernel << " s\n";

    bo_nu1.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu2.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu3.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nu4.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

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

    double t_total = std::chrono::duration<double>(t1 - t_load0).count();
    std::cout << "\n=== Timing breakdown ===\n";
    std::cout << "  CPU extraction + file load : " << t_cpu    << " s\n";
    std::cout << "  k-mer PCIe sync            : " << t_sync   << " s\n";
    std::cout << "  Kernel (HT init + insert)  : " << t_kernel << " s\n";
    std::cout << "  Total end-to-end           : " << t_total  << " s\n";
    std::cout << "  Throughput (vs raw input)  : "
              << (double)(cu_buckets[0].size()+cu_buckets[1].size()+
                          cu_buckets[2].size()+cu_buckets[3].size())
                 * 8.0 / t_kernel / 1e6 << " MB/s (k-mer array)\n";
    return 0;
}
