#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <thread>
#include <stdexcept>

#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

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

/* ---------------------------------------------------------------
 * Routing constants
 *
 * bits[1:0]   → cu_id  (4 CUs per card)
 * bits[9:2]   → seg_id (256 segments per CU)
 * bits[26:10] → T1 index inside kernel (17 bits)   ← NO overlap
 * bits[41:27] → card_id = (h >> CARD_ROUTE_SHIFT) % N_CARD
 *
 * CARD_ROUTE_SHIFT = N_ROUTE_BITS + (SRAM_HT_BITS-1) = 10+17 = 27
 * Using bits above the T1 index avoids clustering within T1.
 * --------------------------------------------------------------- */
static constexpr int N_CARD          = 3;
static constexpr int N_CU_BITS       = 2;
static constexpr int N_SEG_BITS      = 8;
static constexpr int N_ROUTE_BITS    = N_CU_BITS + N_SEG_BITS; /* 10 */
static constexpr int N_CU            = 1 << N_CU_BITS;         /* 4  */
static constexpr int N_SEG           = 1 << N_SEG_BITS;        /* 256 */
static constexpr int N_VIRT_PER_CARD = N_CU * N_SEG;           /* 1024 */
static constexpr int N_VIRT_TOTAL    = N_CARD * N_VIRT_PER_CARD; /* 3072 */
static constexpr int CARD_ROUTE_SHIFT = 27; /* N_ROUTE_BITS + SRAM_HT_BITS(18) - 1 */

static void extract_kmers_fastq(
    const uint8_t*                      seq,
    uint64_t                            len,
    int                                 k,
    uint64_t                            km_mask,
    int                                 shift,
    std::vector<std::vector<uint64_t>>& virt_buckets
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

                    int cu_id   = (int)(h & (N_CU  - 1));
                    int seg_id  = (int)((h >> N_CU_BITS) & (N_SEG - 1));
                    int card_id = (int)((h >> CARD_ROUTE_SHIFT) % (uint64_t)N_CARD);

                    virt_buckets[card_id * N_VIRT_PER_CARD +
                                 cu_id  * N_SEG             +
                                 seg_id].push_back(h);
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

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <xclbin> <fastq1> <fastq2> <k>"
                     " [-d0 <dev>] [-d1 <dev>] [-d2 <dev>]\n";
        return 1;
    }
    std::string xclbin_path = argv[1];
    std::string fq1_path    = argv[2];
    std::string fq2_path    = argv[3];
    int         k           = std::stoi(argv[4]);
    unsigned    dev_idx[N_CARD] = {0, 1, 2};
    for (int i = 5; i < argc - 1; i++) {
        if      (std::string(argv[i]) == "-d0") dev_idx[0] = std::stoul(argv[i+1]);
        else if (std::string(argv[i]) == "-d1") dev_idx[1] = std::stoul(argv[i+1]);
        else if (std::string(argv[i]) == "-d2") dev_idx[2] = std::stoul(argv[i+1]);
    }
    if (k < 1 || k > 31) { std::cerr << "k must be in [1,31]\n"; return 1; }

    uint64_t km_mask = (1ULL << (2*k)) - 1;
    int      shift   = (k-1)*2;

    /* ---- CPU: read FASTQ and route k-mers to 3072 virtual buckets ---- */
    std::cout << "[opt1] Loading FASTQ files ...\n";
    auto t_load0 = std::chrono::high_resolution_clock::now();
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    std::cout << "  FASTQ1: " << fq1_path << "  (" << seq1.size() << " B)\n";
    std::cout << "  FASTQ2: " << fq2_path << "  (" << seq2.size() << " B)\n";

    std::cout << "[opt1] CPU k-mer extraction (k=" << k
              << ", " << N_CARD << " cards × " << N_CU << " CU × "
              << N_SEG << " seg = " << N_VIRT_TOTAL << " virtual buckets) ...\n";

    std::vector<std::vector<uint64_t>> virt_buckets(N_VIRT_TOTAL);
    extract_kmers_fastq(seq1.data(), seq1.size(), k, km_mask, shift, virt_buckets);
    extract_kmers_fastq(seq2.data(), seq2.size(), k, km_mask, shift, virt_buckets);
    seq1.clear(); seq1.shrink_to_fit();
    seq2.clear(); seq2.shrink_to_fit();

    auto t_load1 = std::chrono::high_resolution_clock::now();
    double t_cpu = std::chrono::duration<double>(t_load1 - t_load0).count();

    uint64_t total_kmers = 0;
    for (int c = 0; c < N_CARD; c++) {
        uint64_t card_total = 0;
        for (int cu = 0; cu < N_CU; cu++)
            for (int s = 0; s < N_SEG; s++)
                card_total += virt_buckets[c*N_VIRT_PER_CARD + cu*N_SEG + s].size();
        std::cout << "  Card" << c << " k-mers: " << card_total
                  << " (dev " << dev_idx[c] << ")\n";
        total_kmers += card_total;
    }
    std::cout << "  Total k-mers : " << total_kmers
              << "  (CPU extraction: " << t_cpu << " s)\n";

    /* ---- XRT setup: open 3 devices, load same xclbin ---- */
    std::cout << "[opt1] Opening " << N_CARD << " devices ...\n";
    xrt::device device[N_CARD];
    xrt::uuid   uuid[N_CARD];
    for (int c = 0; c < N_CARD; c++) {
        device[c] = xrt::device(dev_idx[c]);
        uuid[c]   = device[c].load_xclbin(xclbin_path);
        std::cout << "  Card" << c << " (dev " << dev_idx[c] << "): xclbin loaded\n";
    }

    const char* cu_names[N_CU] = {
        "kmer_count_sram_opt1:{kmer_count_sram_opt1_1}",
        "kmer_count_sram_opt1:{kmer_count_sram_opt1_2}",
        "kmer_count_sram_opt1:{kmer_count_sram_opt1_3}",
        "kmer_count_sram_opt1:{kmer_count_sram_opt1_4}"
    };
    xrt::kernel krnl[N_CARD][N_CU];
    for (int c = 0; c < N_CARD; c++)
        for (int cu = 0; cu < N_CU; cu++)
            krnl[c][cu] = xrt::kernel(device[c], uuid[c], cu_names[cu]);

    /* ---- Allocate BOs, build contiguous k-mer arrays ---- */
    std::cout << "[opt1] Allocating BOs ...\n";
    auto t_sync0 = std::chrono::high_resolution_clock::now();

    xrt::bo bo_km[N_CARD][N_CU];
    xrt::bo bo_se[N_CARD][N_CU];
    xrt::bo bo_nu[N_CARD][N_CU];
    uint64_t n_km[N_CARD][N_CU] = {};

    for (int c = 0; c < N_CARD; c++) {
        for (int cu = 0; cu < N_CU; cu++) {
            uint64_t cu_total = 0;
            for (int s = 0; s < N_SEG; s++)
                cu_total += virt_buckets[c*N_VIRT_PER_CARD + cu*N_SEG + s].size();
            n_km[c][cu] = cu_total;

            uint64_t km_bytes = cu_total * sizeof(uint64_t);
            if (km_bytes == 0) km_bytes = 8;

            bo_km[c][cu] = xrt::bo(device[c], km_bytes,      krnl[c][cu].group_id(0));
            bo_se[c][cu] = xrt::bo(device[c], N_SEG*8,       krnl[c][cu].group_id(1));
            bo_nu[c][cu] = xrt::bo(device[c], sizeof(uint64_t), krnl[c][cu].group_id(3));

            uint64_t* km_ptr = bo_km[c][cu].map<uint64_t*>();
            uint64_t* se_ptr = bo_se[c][cu].map<uint64_t*>();
            *bo_nu[c][cu].map<uint64_t*>() = 0;

            uint64_t offset = 0;
            for (int s = 0; s < N_SEG; s++) {
                const auto& bkt = virt_buckets[c*N_VIRT_PER_CARD + cu*N_SEG + s];
                if (!bkt.empty())
                    std::memcpy(km_ptr + offset, bkt.data(), bkt.size() * 8);
                offset += bkt.size();
                se_ptr[s] = offset;
            }
        }
    }
    /* Free routing buckets */
    virt_buckets.clear();
    virt_buckets.shrink_to_fit();

    /* ---- Parallel PCIe sync (one thread per card) ---- */
    std::cout << "[opt1] PCIe sync (parallel across " << N_CARD << " cards) ...\n";
    {
        std::thread sync_t[N_CARD];
        for (int c = 0; c < N_CARD; c++) {
            sync_t[c] = std::thread([&, c]() {
                for (int cu = 0; cu < N_CU; cu++) {
                    bo_km[c][cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);
                    bo_se[c][cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);
                    bo_nu[c][cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);
                }
            });
        }
        for (int c = 0; c < N_CARD; c++) sync_t[c].join();
    }

    auto t_sync1 = std::chrono::high_resolution_clock::now();
    double t_sync = std::chrono::duration<double>(t_sync1 - t_sync0).count();
    std::cout << "  k-mer PCIe sync (parallel): " << t_sync << " s\n";

    /* ---- Launch all 12 CUs simultaneously ---- */
    std::cout << "[opt1] Launching " << N_CARD*N_CU << " CUs ("
              << N_CARD << " cards × " << N_CU << " CU, lazy-clear) ...\n";
    auto t0 = std::chrono::high_resolution_clock::now();

    xrt::run run[N_CARD][N_CU];
    for (int c = 0; c < N_CARD; c++)
        for (int cu = 0; cu < N_CU; cu++)
            run[c][cu] = krnl[c][cu](bo_km[c][cu], bo_se[c][cu],
                                     n_km[c][cu],  bo_nu[c][cu]);
    for (int c = 0; c < N_CARD; c++)
        for (int cu = 0; cu < N_CU; cu++)
            run[c][cu].wait();

    auto t1 = std::chrono::high_resolution_clock::now();
    double t_kernel = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[opt1] All " << N_CARD*N_CU << " CUs finished in " << t_kernel << " s\n";

    /* ---- Collect results ---- */
    uint64_t nu[N_CARD][N_CU] = {};
    for (int c = 0; c < N_CARD; c++) {
        for (int cu = 0; cu < N_CU; cu++) {
            bo_nu[c][cu].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            nu[c][cu] = *bo_nu[c][cu].map<uint64_t*>();
        }
    }
    uint64_t total_unique = 0;
    for (int c = 0; c < N_CARD; c++)
        for (int cu = 0; cu < N_CU; cu++)
            total_unique += nu[c][cu];

    std::cout << "\n=== Results ===\n";
    for (int c = 0; c < N_CARD; c++) {
        uint64_t card_unique = 0;
        for (int cu = 0; cu < N_CU; cu++) card_unique += nu[c][cu];
        std::cout << "  Card" << c << " unique " << k << "-mers: " << card_unique << "\n";
    }
    std::cout << "  Total unique " << k << "-mers: " << total_unique << "\n";

    double t_total = std::chrono::duration<double>(t1 - t_load0).count();
    std::cout << "\n=== Timing breakdown ===\n";
    std::cout << "  CPU extraction + file load      : " << t_cpu    << " s\n";
    std::cout << "  k-mer PCIe sync (parallel)      : " << t_sync   << " s\n";
    std::cout << "  Kernel (" << N_CARD << " cards × " << N_CU
              << " CU, lazy-clear) : " << t_kernel << " s\n";
    std::cout << "  Total end-to-end                : " << t_total  << " s\n";
    std::cout << "  Throughput (k-mer array)        : "
              << (double)total_kmers * 8.0 / t_kernel / 1e6 << " MB/s\n";
    return 0;
}
