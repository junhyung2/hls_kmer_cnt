#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <stdexcept>

#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

/* ---------------------------------------------------------------
 * NT4 lookup table and canonical k-mer extraction
 * (same as opt1/card2)
 * --------------------------------------------------------------- */
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
 * Routing constants (must match kernel header)
 * --------------------------------------------------------------- */
static constexpr int N_CU_BITS  = 2;
static constexpr int N_SEG_BITS = 8;              /* must match kmer_count_sram.h */
static constexpr int N_ROUTE_BITS = N_CU_BITS + N_SEG_BITS;  /* 10 */
static constexpr int N_CU  = 1 << N_CU_BITS;   /* 4   */
static constexpr int N_SEG = 1 << N_SEG_BITS;  /* 256 */
static constexpr int N_VIRT = N_CU * N_SEG;    /* 1024 virtual buckets */

/* ---------------------------------------------------------------
 * Extract canonical k-mer hashes and route to 64 virtual buckets:
 *   cu_id  = h & (N_CU  - 1)             [bits 0-1]
 *   seg_id = (h >> N_CU_BITS) & (N_SEG - 1)  [bits 2-5]
 *   virt   = cu_id | (seg_id << N_CU_BITS)
 * --------------------------------------------------------------- */
static void extract_kmers_fastq(
    const uint8_t*                         seq,
    uint64_t                               len,
    int                                    k,
    uint64_t                               km_mask,
    int                                    shift,
    std::vector<std::vector<uint64_t>>&    virt_buckets  /* [N_VIRT=64] */
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
                    int cu_id  = (int)(h & (N_CU  - 1));
                    int seg_id = (int)((h >> N_CU_BITS) & (N_SEG - 1));
                    virt_buckets[cu_id * N_SEG + seg_id].push_back(h);
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

/* ---------------------------------------------------------------
 * HBM host-side size for seg_ends and n_unique BOs (must match
 * kernel group_id mapping in connectivity cfg).
 * We pass N_SEG uint64_t values as seg_ends → 128 bytes.
 * --------------------------------------------------------------- */
#ifndef SRAM_HT_BITS_HOST
  #define SRAM_HT_BITS_HOST 18
#endif

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <xclbin> <fastq1> <fastq2> <k> [-d <dev_idx>]\n";
        return 1;
    }
    std::string xclbin_path = argv[1];
    std::string fq1_path    = argv[2];
    std::string fq2_path    = argv[3];
    int         k           = std::stoi(argv[4]);
    unsigned    dev_idx     = 0;
    for (int i = 5; i < argc - 1; i++)
        if (std::string(argv[i]) == "-d") dev_idx = std::stoul(argv[i+1]);
    if (k < 1 || k > 31) { std::cerr << "k must be in [1,31]\n"; return 1; }

    uint64_t km_mask = (1ULL << (2*k)) - 1;
    int      shift   = (k-1)*2;

    /* ---- CPU: read FASTQ and route k-mers to 64 virtual buckets ---- */
    std::cout << "[sram] Loading FASTQ files ...\n";
    auto t_load0 = std::chrono::high_resolution_clock::now();
    auto seq1 = read_file(fq1_path);
    auto seq2 = read_file(fq2_path);
    std::cout << "  FASTQ1: " << fq1_path << "  (" << seq1.size() << " B)\n";
    std::cout << "  FASTQ2: " << fq2_path << "  (" << seq2.size() << " B)\n";

    std::cout << "[sram] CPU k-mer extraction (k=" << k
              << ", " << N_VIRT << " virtual buckets) ...\n";
    /* virt_buckets[cu_id * N_SEG + seg_id] */
    std::vector<std::vector<uint64_t>> virt_buckets(N_VIRT);
    extract_kmers_fastq(seq1.data(), seq1.size(), k, km_mask, shift, virt_buckets);
    extract_kmers_fastq(seq2.data(), seq2.size(), k, km_mask, shift, virt_buckets);
    seq1.clear(); seq1.shrink_to_fit();
    seq2.clear(); seq2.shrink_to_fit();

    auto t_load1 = std::chrono::high_resolution_clock::now();
    double t_cpu = std::chrono::duration<double>(t_load1 - t_load0).count();

    uint64_t total_kmers = 0;
    for (int cu = 0; cu < N_CU; cu++) {
        uint64_t cu_total = 0;
        for (int s = 0; s < N_SEG; s++) cu_total += virt_buckets[cu * N_SEG + s].size();
        std::cout << "  CU" << cu << " k-mers: " << cu_total
                  << " (" << N_SEG << " segments)\n";
        total_kmers += cu_total;
    }
    std::cout << "  Total k-mers : " << total_kmers
              << "  (CPU extraction: " << t_cpu << " s)\n";

    /* ---- XRT setup ---- */
    std::cout << "[sram] Opening device " << dev_idx << " ...\n";
    xrt::device device(dev_idx);
    xrt::uuid   uuid = device.load_xclbin(xclbin_path);

    const char* cu_names[4] = {
        "kmer_count_sram:{kmer_count_sram_1}",
        "kmer_count_sram:{kmer_count_sram_2}",
        "kmer_count_sram:{kmer_count_sram_3}",
        "kmer_count_sram:{kmer_count_sram_4}"
    };
    xrt::kernel krnl[N_CU];
    for (int cu = 0; cu < N_CU; cu++)
        krnl[cu] = xrt::kernel(device, uuid, cu_names[cu]);

    /* ---- Allocate BOs and sync k-mers ---- */
    std::cout << "[sram] Allocating BOs and syncing k-mers ...\n";
    auto t_sync0 = std::chrono::high_resolution_clock::now();

    xrt::bo bo_km[N_CU];
    xrt::bo bo_se[N_CU];   /* seg_ends: N_SEG × 8 bytes */
    xrt::bo bo_nu[N_CU];   /* n_unique: 8 bytes          */

    for (int cu = 0; cu < N_CU; cu++) {
        /* Build contiguous k-mer array (seg0 | seg1 | ... | seg15) */
        uint64_t cu_total = 0;
        for (int s = 0; s < N_SEG; s++) cu_total += virt_buckets[cu * N_SEG + s].size();
        uint64_t km_bytes = cu_total * sizeof(uint64_t);
        if (km_bytes == 0) km_bytes = 8;

        bo_km[cu] = xrt::bo(device, km_bytes, krnl[cu].group_id(0));
        uint64_t* km_ptr = bo_km[cu].map<uint64_t*>();

        /* seg_ends BO: N_SEG uint64_t values */
        bo_se[cu] = xrt::bo(device, N_SEG * sizeof(uint64_t), krnl[cu].group_id(1));
        uint64_t* se_ptr = bo_se[cu].map<uint64_t*>();

        uint64_t offset = 0;
        for (int s = 0; s < N_SEG; s++) {
            const auto& bkt = virt_buckets[cu * N_SEG + s];
            if (!bkt.empty())
                std::memcpy(km_ptr + offset, bkt.data(), bkt.size() * sizeof(uint64_t));
            offset += bkt.size();
            se_ptr[s] = offset;          /* cumulative end offset */
        }

        bo_km[cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);
        bo_se[cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);

        bo_nu[cu] = xrt::bo(device, sizeof(uint64_t), krnl[cu].group_id(3));
        *bo_nu[cu].map<uint64_t*>() = 0;
        bo_nu[cu].sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }

    auto t_sync1 = std::chrono::high_resolution_clock::now();
    double t_sync = std::chrono::duration<double>(t_sync1 - t_sync0).count();
    std::cout << "  k-mer PCIe sync: " << t_sync << " s\n";

    /* ---- Launch all 4 CUs simultaneously ---- */
    std::cout << "[sram] Launching " << N_CU << " CUs (URAM multi-segment) ...\n";
    auto t0 = std::chrono::high_resolution_clock::now();

    /* Read n_kmers per CU from seg_ends[N_SEG-1] (already in host-mapped BO) */
    uint64_t n_km[N_CU];
    for (int cu = 0; cu < N_CU; cu++)
        n_km[cu] = bo_se[cu].map<uint64_t*>()[N_SEG - 1];

    xrt::run run[N_CU];
    for (int cu = 0; cu < N_CU; cu++)
        run[cu] = krnl[cu](bo_km[cu], bo_se[cu], n_km[cu], bo_nu[cu]);
    for (int cu = 0; cu < N_CU; cu++) run[cu].wait();

    auto t1 = std::chrono::high_resolution_clock::now();
    double t_kernel = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[sram] All " << N_CU << " CUs finished in " << t_kernel << " s\n";

    /* ---- Collect results ---- */
    uint64_t nu[N_CU] = {};
    for (int cu = 0; cu < N_CU; cu++) {
        bo_nu[cu].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        nu[cu] = *bo_nu[cu].map<uint64_t*>();
    }
    uint64_t total_unique = 0;
    for (int cu = 0; cu < N_CU; cu++) total_unique += nu[cu];

    std::cout << "\n=== Results ===\n";
    for (int cu = 0; cu < N_CU; cu++)
        std::cout << "  CU" << cu << " unique " << k << "-mers: " << nu[cu] << "\n";
    std::cout << "  Total unique " << k << "-mers: " << total_unique << "\n";

    double t_total = std::chrono::duration<double>(t1 - t_load0).count();
    std::cout << "\n=== Timing breakdown ===\n";
    std::cout << "  CPU extraction + file load : " << t_cpu    << " s\n";
    std::cout << "  k-mer PCIe sync            : " << t_sync   << " s\n";
    std::cout << "  Kernel (URAM multi-seg)    : " << t_kernel << " s\n";
    std::cout << "  Total end-to-end           : " << t_total  << " s\n";
    std::cout << "  Throughput (k-mer array)   : "
              << (double)total_kmers * 8.0 / t_kernel / 1e6 << " MB/s\n";
    return 0;
}
