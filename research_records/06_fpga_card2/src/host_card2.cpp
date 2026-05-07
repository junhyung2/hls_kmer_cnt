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

/* Extract canonical k-mer hashes, routing to 8 buckets (h & 7) for 2 cards × 4 CUs */
static void extract_kmers_fastq(
    const uint8_t*                      seq,
    uint64_t                            len,
    int                                 k,
    uint64_t                            km_mask,
    int                                 shift,
    std::vector<std::vector<uint64_t>>& cu_buckets  /* [8] */
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
                    cu_buckets[h & 7].push_back(h);
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
static const int   N_CU_PER_CARD = 4;
static const int   N_CARDS       = 2;
static const int   N_CU_TOTAL    = N_CU_PER_CARD * N_CARDS;

struct DeviceContext {
    xrt::device device;
    xrt::uuid   uuid;
    xrt::kernel krnl[4];
    xrt::bo     bo_km[4];
    xrt::bo     bo_t1[4];
    xrt::bo     bo_t2[4];
    xrt::bo     bo_nu[4];
    xrt::run    run[4];
};

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <xclbin> <fastq1> <fastq2> <k> [-d0 <dev0>] [-d1 <dev1>]\n";
        return 1;
    }
    std::string xclbin_path = argv[1];
    std::string fq1_path    = argv[2];
    std::string fq2_path    = argv[3];
    int         k           = std::stoi(argv[4]);
    unsigned    dev0_idx    = 0;
    unsigned    dev1_idx    = 1;
    for (int i = 5; i < argc - 1; i++) {
        if (std::string(argv[i]) == "-d0") dev0_idx = std::stoul(argv[i+1]);
        if (std::string(argv[i]) == "-d1") dev1_idx = std::stoul(argv[i+1]);
    }
    if (k < 1 || k > 31) { std::cerr << "k must be in [1,31]\n"; return 1; }

    uint64_t km_mask = (1ULL << (2*k)) - 1;
    int      shift   = (k-1)*2;

    /* ---- CPU: read FASTQ and extract k-mers to 8 buckets ---- */
    std::cout << "[card2] Loading FASTQ files ...\n";
    auto t_load0 = std::chrono::high_resolution_clock::now();
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    std::cout << "  FASTQ1: " << fq1_path << "  (" << seq1.size() << " B)\n";
    std::cout << "  FASTQ2: " << fq2_path << "  (" << seq2.size() << " B)\n";

    std::cout << "[card2] CPU k-mer extraction (k=" << k << ", 8 buckets) ...\n";
    std::vector<std::vector<uint64_t>> cu_buckets(N_CU_TOTAL);
    extract_kmers_fastq(seq1.data(), seq1.size(), k, km_mask, shift, cu_buckets);
    extract_kmers_fastq(seq2.data(), seq2.size(), k, km_mask, shift, cu_buckets);
    seq1.clear(); seq1.shrink_to_fit();
    seq2.clear(); seq2.shrink_to_fit();

    auto t_load1 = std::chrono::high_resolution_clock::now();
    double t_cpu = std::chrono::duration<double>(t_load1 - t_load0).count();
    uint64_t total_kmers = 0;
    for (int i = 0; i < N_CU_TOTAL; i++) {
        std::cout << "  CU" << i << " (card" << (i/4) << ") k-mers: " << cu_buckets[i].size() << "\n";
        total_kmers += cu_buckets[i].size();
    }
    std::cout << "  Total k-mers: " << total_kmers
              << "  (CPU extraction: " << t_cpu << " s)\n";

    /* ---- XRT setup: open both devices ---- */
    std::cout << "[card2] Opening device " << dev0_idx << " (card0) and "
              << dev1_idx << " (card1) ...\n";

    DeviceContext ctx[2];
    unsigned dev_idx[2] = { dev0_idx, dev1_idx };
    const char* cu_names[4] = {
        "kmer_count_opt1:{kmer_count_opt1_1}",
        "kmer_count_opt1:{kmer_count_opt1_2}",
        "kmer_count_opt1:{kmer_count_opt1_3}",
        "kmer_count_opt1:{kmer_count_opt1_4}"
    };

    for (int c = 0; c < N_CARDS; c++) {
        ctx[c].device = xrt::device(dev_idx[c]);
        ctx[c].uuid   = ctx[c].device.load_xclbin(xclbin_path);
        for (int cu = 0; cu < N_CU_PER_CARD; cu++)
            ctx[c].krnl[cu] = xrt::kernel(ctx[c].device, ctx[c].uuid, cu_names[cu]);
    }

    /* ---- Allocate + sync k-mer buffers, allocate HT buffers ---- */
    std::cout << "[card2] Allocating k-mer buffers and syncing ...\n";
    auto t_sync0 = std::chrono::high_resolution_clock::now();

    for (int c = 0; c < N_CARDS; c++) {
        for (int cu = 0; cu < N_CU_PER_CARD; cu++) {
            int bucket = c * N_CU_PER_CARD + cu;
            uint64_t nbytes = cu_buckets[bucket].size() * sizeof(uint64_t);
            if (nbytes == 0) nbytes = 8;
            ctx[c].bo_km[cu] = xrt::bo(ctx[c].device, nbytes, ctx[c].krnl[cu].group_id(0));
            if (!cu_buckets[bucket].empty())
                std::memcpy(ctx[c].bo_km[cu].map<uint64_t*>(),
                            cu_buckets[bucket].data(), nbytes);
            ctx[c].bo_km[cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);

            ctx[c].bo_t1[cu] = xrt::bo(ctx[c].device, HT_T_BYTES, ctx[c].krnl[cu].group_id(1));
            ctx[c].bo_t2[cu] = xrt::bo(ctx[c].device, HT_T_BYTES, ctx[c].krnl[cu].group_id(2));

            ctx[c].bo_nu[cu] = xrt::bo(ctx[c].device, sizeof(uint64_t), ctx[c].krnl[cu].group_id(4));
            *ctx[c].bo_nu[cu].map<uint64_t*>() = 0;
            ctx[c].bo_nu[cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);
        }
    }

    auto t_sync1 = std::chrono::high_resolution_clock::now();
    double t_sync = std::chrono::duration<double>(t_sync1 - t_sync0).count();
    std::cout << "  k-mer sync (both cards): " << t_sync << " s\n";

    /* ---- Launch all 8 CUs across 2 cards simultaneously ---- */
    std::cout << "[card2] Launching 8 CUs (4 per card, both cards in parallel) ...\n";
    auto t0 = std::chrono::high_resolution_clock::now();

    for (int c = 0; c < N_CARDS; c++) {
        for (int cu = 0; cu < N_CU_PER_CARD; cu++) {
            int bucket = c * N_CU_PER_CARD + cu;
            ctx[c].run[cu] = ctx[c].krnl[cu](
                ctx[c].bo_km[cu],
                ctx[c].bo_t1[cu],
                ctx[c].bo_t2[cu],
                (uint64_t)cu_buckets[bucket].size(),
                ctx[c].bo_nu[cu]
            );
        }
    }

    /* Wait for all CUs on both cards */
    for (int c = 0; c < N_CARDS; c++)
        for (int cu = 0; cu < N_CU_PER_CARD; cu++)
            ctx[c].run[cu].wait();

    auto t1 = std::chrono::high_resolution_clock::now();
    double t_kernel = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[card2] All kernels finished in " << t_kernel << " s\n";

    /* ---- Collect results ---- */
    uint64_t nu[N_CU_TOTAL] = {};
    for (int c = 0; c < N_CARDS; c++) {
        for (int cu = 0; cu < N_CU_PER_CARD; cu++) {
            ctx[c].bo_nu[cu].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            nu[c * N_CU_PER_CARD + cu] = *ctx[c].bo_nu[cu].map<uint64_t*>();
        }
    }

    uint64_t total_unique = 0;
    for (int i = 0; i < N_CU_TOTAL; i++) total_unique += nu[i];

    std::cout << "\n=== Results ===\n";
    for (int i = 0; i < N_CU_TOTAL; i++)
        std::cout << "  CU" << i << " (card" << (i/4) << ") unique " << k << "-mers : " << nu[i] << "\n";
    std::cout << "  Total unique " << k << "-mers : " << total_unique << "\n";

    double t_total = std::chrono::duration<double>(t1 - t_load0).count();
    std::cout << "\n=== Timing breakdown ===\n";
    std::cout << "  CPU extraction + file load : " << t_cpu    << " s\n";
    std::cout << "  k-mer PCIe sync (2 cards)  : " << t_sync   << " s\n";
    std::cout << "  Kernel (both cards)        : " << t_kernel << " s\n";
    std::cout << "  Total end-to-end           : " << t_total  << " s\n";
    std::cout << "  Throughput (k-mer array)   : "
              << (double)total_kmers * 8.0 / t_kernel / 1e6 << " MB/s\n";
    return 0;
}
