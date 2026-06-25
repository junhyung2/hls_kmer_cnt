/* exp58: host_3card_sort58.cpp — FPGA extraction (no hash) + sub-bin scatter sort
 *
 * vs exp55 (radix sort):
 *   exp55: kmer_pool[N_CU] flat array → sort entire pool (144GB for 50G)
 *          → 8-CU parallel impossible (2× pool would need 288GB RAM)
 *          → sequential sort: 185s wall
 *
 *   exp58: km_bins[N_CU][N_SBINS] — scatter to 64 sub-bins during append
 *          → 8-CU always parallel (peak RAM = pool + 8 × max_sub_bin temp)
 *          → per-sub-bin sort: ~280MB each, better cache utilization
 *          FPGA: outputs canonical k-mer (no hash), 1-XOR CU routing
 *
 * Sub-bin assignment: sbin = (canonical >> SB_SHIFT) & SB_MASK
 *   SB_SHIFT=3, SB_MASK=0x3F → bits [3:8] of the 42-bit canonical value
 *   Deduplication: same canonical → same (CU, sbin) guaranteed ✓
 *
 * Memory (50G): pool=144GB (sub-bin arrays) + sort_temp=8×280MB + flat_buf=2.3GB ≈ 149GB
 */

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <thread>
#include <stdexcept>
#include <array>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

#include "fastq_extractor_8port.h"

static constexpr int N_DEV          = 3;
static constexpr int N_INST_PER_DEV = 3;
static constexpr int N_INST_TOTAL   = N_DEV * N_INST_PER_DEV;
static constexpr int N_CU           = FE_N_CU;

static constexpr int      KM_BITS   = 42;
static constexpr uint64_t KEY_MASK  = (1ULL << KM_BITS) - 1;

/* Sub-bin parameters: 64 sub-bins per CU, using bits [3:8] of canonical value */
static constexpr int      N_SBINS   = 64;
static constexpr int      SB_BITS   = 6;
static constexpr int      SB_SHIFT  = 3;
static constexpr uint64_t SB_MASK   = (1ULL << SB_BITS) - 1;

struct MmapFile {
    const uint8_t* data = nullptr;
    uint64_t       size = 0;
    int            fd   = -1;

    explicit MmapFile(const std::string& path) {
        fd = open(path.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("Cannot open: " + path);
        struct stat st;
        if (fstat(fd, &st) < 0) throw std::runtime_error("fstat failed: " + path);
        size = (uint64_t)st.st_size;
        if (size > 0) {
            data = static_cast<const uint8_t*>(
                mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
            if (data == MAP_FAILED) {
                data = nullptr;
                throw std::runtime_error("mmap failed: " + path);
            }
            madvise((void*)data, size, MADV_SEQUENTIAL);
        }
    }
    ~MmapFile() {
        if (data) munmap((void*)data, size);
        if (fd >= 0) close(fd);
    }
    MmapFile(const MmapFile&) = delete;
    MmapFile& operator=(const MmapFile&) = delete;
};

static inline uint8_t cat_byte(const MmapFile& f1, const MmapFile& f2, uint64_t pos)
{
    return pos < f1.size ? f1.data[pos] : f2.data[pos - f1.size];
}

static uint64_t find_next_header_at(const MmapFile& f1, const MmapFile& f2,
                                     uint64_t total, uint64_t target)
{
    uint64_t i = (target > 0) ? target - 1 : 0;
    for (; i + 1 < total; i++) {
        if (cat_byte(f1, f2, i) != '\n') continue;
        if (cat_byte(f1, f2, i + 1) != '@') continue;
        uint64_t p = i + 2;
        while (p < total && cat_byte(f1, f2, p) != '\n') p++;
        if (p >= total) break;
        p++;
        int bases = 0;
        bool valid = true;
        while (p < total && cat_byte(f1, f2, p) != '\n') {
            uint8_t c = cat_byte(f1, f2, p++);
            if (c=='A'||c=='C'||c=='G'||c=='T'||c=='N'||
                c=='a'||c=='c'||c=='g'||c=='t'||c=='n') { bases++; }
            else { valid = false; break; }
        }
        if (valid && bases >= 10) return i + 1;
    }
    return total;
}

static void cat_memcpy(uint8_t* dst, const MmapFile& f1, const MmapFile& f2,
                        uint64_t off, uint64_t len)
{
    if (off + len <= f1.size) {
        std::memcpy(dst, f1.data + off, len);
    } else if (off >= f1.size) {
        std::memcpy(dst, f2.data + (off - f1.size), len);
    } else {
        uint64_t part1 = f1.size - off;
        std::memcpy(dst, f1.data + off, part1);
        std::memcpy(dst + part1, f2.data, len - part1);
    }
}

static double elapsed(std::chrono::high_resolution_clock::time_point t0,
                      std::chrono::high_resolution_clock::time_point t1)
{
    return std::chrono::duration<double>(t1 - t0).count();
}

static inline uint64_t hash64_host(uint64_t key, uint64_t mask)
{
    key ^= key >> 33;
    key ^= key << 21;
    return key & mask;
}

/* 3-pass LSD radix sort for 42-bit keys.
 * After 3 passes (odd), sorted result is in dst[].
 * src[] is overwritten (used as ping-pong buffer) — caller must not reuse src after call.
 */
static void lsd_radix_sort_42(const uint64_t* src, uint64_t* dst, uint64_t n)
{
    constexpr int      RB    = 14;
    constexpr uint64_t RADIX = 1ULL << RB;
    constexpr uint64_t MASK  = RADIX - 1;
    uint64_t* rd = const_cast<uint64_t*>(src);
    uint64_t* wr = dst;
    for (int pass = 0; pass < 3; pass++) {
        const int shift = pass * RB;
        uint64_t cnt[RADIX] = {};
        for (uint64_t i = 0; i < n; i++) cnt[(rd[i] >> shift) & MASK]++;
        uint64_t psum[RADIX];
        psum[0] = 0;
        for (uint64_t i = 1; i < RADIX; i++) psum[i] = psum[i-1] + cnt[i-1];
        for (uint64_t i = 0; i < n; i++) wr[psum[(rd[i] >> shift) & MASK]++] = rd[i];
        std::swap(rd, wr);
    }
    /* After 3 passes: sorted result is in dst (rd points to dst after last swap) */
}

static uint64_t count_sorted(const uint64_t* arr, uint64_t n, uint64_t* histogram)
{
    if (n == 0) return 0;
    uint64_t unique = 0, prev = arr[0], cnt = 1;
    for (uint64_t i = 1; i < n; i++) {
        if (arr[i] == prev) { cnt++; }
        else {
            histogram[cnt < 255 ? cnt : 255]++;
            unique++;
            prev = arr[i];
            cnt  = 1;
        }
    }
    histogram[cnt < 255 ? cnt : 255]++;
    return unique + 1;
}

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <xclbin> <fastq1> <fastq2> <k>"
                     " [-b <mb_per_inst>] [-n <runs>] [-hist]\n";
        return 1;
    }
    std::string xclbin_path = argv[1];
    std::string fq1_path    = argv[2];
    std::string fq2_path    = argv[3];
    int         k           = std::stoi(argv[4]);
    uint64_t    batch_mb    = 64;
    int         n_runs      = 1;
    bool        print_hist  = false;
    for (int i = 5; i < argc; i++) {
        if (std::string(argv[i]) == "-b" && i+1<argc) batch_mb = std::stoull(argv[++i]);
        if (std::string(argv[i]) == "-n" && i+1<argc) n_runs   = std::stoi(argv[++i]);
        if (std::string(argv[i]) == "-hist") print_hist = true;
    }
    if (k != FE_K) { std::cerr << "Only k=" << FE_K << " supported\n"; return 1; }

    MmapFile fq1(fq1_path), fq2(fq2_path);
    uint64_t fq_total = fq1.size + fq2.size;
    std::cout << "[sort58] Total FASTQ: " << fq_total / 1048576.0 << " MB\n";

    uint64_t batch_bytes = batch_mb * 1048576ULL;
    int est_batches = (int)((fq_total + (uint64_t)N_INST_TOTAL * batch_bytes - 1)
                            / ((uint64_t)N_INST_TOTAL * batch_bytes));
    std::cout << "  batch_bytes/inst: " << batch_bytes / 1048576.0 << " MB"
              << "  estimated batches: " << est_batches << "\n";

    /* ── Open devices ── */
    xrt::device device[N_DEV];
    xrt::uuid   uuid[N_DEV];
    for (int d = 0; d < N_DEV; d++) {
        device[d] = xrt::device(d);
        uuid[d]   = device[d].load_xclbin(xclbin_path);
    }

    const char* krnl_names[N_INST_PER_DEV] = {
        "fastq_extractor_8port:{fastq_extractor_8port_1}",
        "fastq_extractor_8port:{fastq_extractor_8port_2}",
        "fastq_extractor_8port:{fastq_extractor_8port_3}",
    };

    /* ── Allocate BOs ── */
    const uint64_t cu_bytes  = (uint64_t)FE_MAX_PER_CU * sizeof(uint64_t);
    const uint64_t se_bytes  = (uint64_t)N_CU * sizeof(uint64_t);
    const uint64_t fq_bytes  = (batch_bytes + 8192 + 7) & ~7ULL;

    xrt::kernel ext_krnl[N_DEV][N_INST_PER_DEV];
    xrt::bo     fq_bo   [N_DEV][N_INST_PER_DEV][2];
    xrt::bo     km_bo   [N_DEV][N_INST_PER_DEV][N_CU];
    xrt::bo     se_bo   [N_DEV][N_INST_PER_DEV];

    for (int d = 0; d < N_DEV; d++) {
        for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
            ext_krnl[d][ci] = xrt::kernel(device[d], uuid[d], krnl_names[ci]);
            for (int pp = 0; pp < 2; pp++)
                fq_bo[d][ci][pp] = xrt::bo(device[d], fq_bytes,
                                            ext_krnl[d][ci].group_id(0));
            se_bo[d][ci] = xrt::bo(device[d], se_bytes,
                                    ext_krnl[d][ci].group_id(10));
            for (int cu = 0; cu < N_CU; cu++)
                km_bo[d][ci][cu] = xrt::bo(device[d], cu_bytes,
                                            ext_krnl[d][ci].group_id(2 + cu));
        }
    }

    /* ── Flat DMA staging buffers (per CU, interleaved by instance) ── */
    const uint64_t flat_per_cu = (uint64_t)FE_MAX_PER_CU * N_INST_TOTAL;
    std::vector<uint64_t> flat_buf[N_CU];
    for (int cu = 0; cu < N_CU; cu++)
        flat_buf[cu].resize(flat_per_cu, 0);
    uint64_t flat_cnt[N_INST_TOTAL][N_CU] = {};

    /* ── Sub-bin k-mer storage: km_bins[cu][sbin] ── */
    std::vector<uint64_t> km_bins[N_CU][N_SBINS];
    {
        double   est_reads    = fq_total / 300.0;
        uint64_t est_per_cu   = (uint64_t)(est_reads * (150 - k + 1) / N_CU * 1.2);
        uint64_t est_per_sbin = est_per_cu / N_SBINS;
        for (int cu = 0; cu < N_CU; cu++)
            for (int s = 0; s < N_SBINS; s++)
                km_bins[cu][s].reserve(est_per_sbin);
        std::cout << "  Pool reserved: ~" << est_per_cu << " slots/CU ("
                  << est_per_cu * 8.0 / 1e9 << " GB/CU), "
                  << est_per_sbin << " slots/sub-bin\n";
    }

    /* ── Batch offset helpers (same as exp55) ── */
    struct BatchInfo {
        uint64_t start, end, size;
        uint64_t sub_start[N_INST_TOTAL], sub_len[N_INST_TOTAL];
    };

    auto compute_batch_info = [&](uint64_t from) -> BatchInfo {
        BatchInfo bi;
        bi.start = from;
        bi.end   = std::min(from + (uint64_t)N_INST_TOTAL * batch_bytes, fq_total);
        bi.size  = bi.end - bi.start;
        uint64_t cursor = bi.start;
        for (int inst = 0; inst < N_INST_TOTAL; inst++) {
            uint64_t target = bi.start + (uint64_t)(inst + 1) * batch_bytes;
            target = std::min(target, bi.end);
            uint64_t end_i  = (target < fq_total)
                              ? find_next_header_at(fq1, fq2, fq_total, target)
                              : bi.end;
            bi.sub_start[inst] = cursor;
            bi.sub_len[inst]   = end_i - cursor;
            cursor             = end_i;
        }
        bi.size = cursor - bi.start;
        bi.end  = cursor;
        return bi;
    };

    auto fill_and_sync_pp = [&](const BatchInfo& bi, int pp) {
        for (int inst = 0; inst < N_INST_TOTAL; inst++) {
            int d  = inst / N_INST_PER_DEV;
            int ci = inst % N_INST_PER_DEV;
            uint8_t* dst = fq_bo[d][ci][pp].map<uint8_t*>();
            cat_memcpy(dst, fq1, fq2, bi.sub_start[inst], bi.sub_len[inst]);
            fq_bo[d][ci][pp].sync(XCL_BO_SYNC_BO_TO_DEVICE,
                                   bi.sub_len[inst], 0);
        }
    };

    std::vector<double> run_e2e_v;

    for (int run = 0; run < n_runs; run++) {
        std::cout << "\n[sort58] ====== Run " << (run + 1) << " / " << n_runs << " ======\n";

        for (int cu = 0; cu < N_CU; cu++)
            for (int s = 0; s < N_SBINS; s++)
                km_bins[cu][s].clear();

        auto t_run_start = std::chrono::high_resolution_clock::now();
        double sum_ka = 0, sum_append = 0;
        uint64_t total_fastq_processed = 0;
        int      batch_num = 0;

        BatchInfo cur_bi = compute_batch_info(0);
        fill_and_sync_pp(cur_bi, 0);
        int pp_cur = 0;

        while (cur_bi.start < fq_total) {
            batch_num++;
            total_fastq_processed += cur_bi.size;
            std::fill(&flat_cnt[0][0], &flat_cnt[0][0] + N_INST_TOTAL * N_CU, 0ULL);

            auto t_b0 = std::chrono::high_resolution_clock::now();

            int pp_nxt = 1 - pp_cur;
            BatchInfo nxt_bi;
            bool has_next = (cur_bi.end < fq_total);
            std::thread pf_th;
            if (has_next) {
                pf_th = std::thread([&, pp_nxt]() {
                    nxt_bi = compute_batch_info(cur_bi.end);
                    fill_and_sync_pp(nxt_bi, pp_nxt);
                });
            }

            /* Instance workers: launch FPGA, wait, DMA back into flat_buf */
            std::vector<std::thread> inst_workers;
            xrt::run ext_run[N_DEV][N_INST_PER_DEV];
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
                    int inst = d * N_INST_PER_DEV + ci;
                    ext_run[d][ci] = ext_krnl[d][ci](
                        fq_bo[d][ci][pp_cur], cur_bi.sub_len[inst],
                        km_bo[d][ci][0], km_bo[d][ci][1],
                        km_bo[d][ci][2], km_bo[d][ci][3],
                        km_bo[d][ci][4], km_bo[d][ci][5],
                        km_bo[d][ci][6], km_bo[d][ci][7],
                        se_bo[d][ci]);
                }
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
                    inst_workers.emplace_back([&, d, ci]() {
                        int inst = d * N_INST_PER_DEV + ci;
                        ext_run[d][ci].wait();
                        se_bo[d][ci].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
                        const uint64_t* se = se_bo[d][ci].map<const uint64_t*>();
                        for (int cu = 0; cu < N_CU; cu++)
                            flat_cnt[inst][cu] = se[cu];
                        for (int cu = 0; cu < N_CU; cu++) {
                            uint64_t n = flat_cnt[inst][cu];
                            if (n > 0) {
                                km_bo[d][ci][cu].sync(XCL_BO_SYNC_BO_FROM_DEVICE,
                                                       n * sizeof(uint64_t), 0);
                                const uint64_t* src = km_bo[d][ci][cu].map<const uint64_t*>();
                                uint64_t dst_off = (uint64_t)inst * FE_MAX_PER_CU;
                                std::memcpy(flat_buf[cu].data() + dst_off,
                                            src, n * sizeof(uint64_t));
                            }
                        }
                    });
                }

            for (auto& t : inst_workers) t.join();
            auto t_ka1 = std::chrono::high_resolution_clock::now();
            sum_ka += elapsed(t_b0, t_ka1);

            if (pf_th.joinable()) pf_th.join();

            /* Scatter: 8-CU parallel, runs after all DMA is done */
            std::vector<std::thread> scatter_threads;
            for (int cu = 0; cu < N_CU; cu++) {
                scatter_threads.emplace_back([&, cu]() {
                    uint64_t write_off = 0;
                    for (int inst = 0; inst < N_INST_TOTAL; inst++) {
                        uint64_t n = flat_cnt[inst][cu];
                        if (n == 0) continue;
                        uint64_t src_off = (uint64_t)inst * FE_MAX_PER_CU;
                        if (write_off != src_off)
                            std::memmove(flat_buf[cu].data() + write_off,
                                         flat_buf[cu].data() + src_off,
                                         n * sizeof(uint64_t));
                        write_off += n;
                    }
                    const uint64_t* src = flat_buf[cu].data();
                    for (uint64_t i = 0; i < write_off; i++) {
                        uint64_t y    = src[i] & KEY_MASK;
                        uint64_t h    = hash64_host(y, KEY_MASK);
                        int      sbin = (int)((h >> SB_SHIFT) & SB_MASK);
                        km_bins[cu][sbin].push_back(y);
                    }
                });
            }
            for (auto& t : scatter_threads) t.join();
            auto t_app1 = std::chrono::high_resolution_clock::now();
            sum_append += elapsed(t_ka1, t_app1);

            uint64_t pool_total = 0;
            for (int cu = 0; cu < N_CU; cu++)
                for (int s = 0; s < N_SBINS; s++)
                    pool_total += km_bins[cu][s].size();

            std::cout << "  [batch " << batch_num << "]"
                      << "  FASTQ: " << cur_bi.size / 1048576.0 << " MB"
                      << "  KA+DMA: " << elapsed(t_b0, t_ka1) << "s"
                      << "  append: " << elapsed(t_ka1, t_app1) << "s"
                      << "  pool: " << pool_total << " k-mers\n";
            std::cout.flush();

            if (!has_next) break;
            pp_cur = pp_nxt;
            cur_bi = nxt_bi;
        }

        auto t_extract_end = std::chrono::high_resolution_clock::now();
        double t_extract = elapsed(t_run_start, t_extract_end);

        /* ── Sort + count: 8 CU parallel, 64 sub-bins sequential per CU ── */
        uint64_t total_pool = 0;
        uint64_t max_sbin_global = 0;
        for (int cu = 0; cu < N_CU; cu++)
            for (int s = 0; s < N_SBINS; s++) {
                uint64_t sz = km_bins[cu][s].size();
                total_pool      += sz;
                if (sz > max_sbin_global) max_sbin_global = sz;
            }

        std::cout << "\n  Total k-mers in pool: " << total_pool
                  << " (" << total_pool * 8.0 / 1e9 << " GB)\n"
                  << "  Max sub-bin size: " << max_sbin_global
                  << " (" << max_sbin_global * 8.0 / 1e6 << " MB)\n"
                  << "  Sort mode: parallel (8 CU × 64 sub-bins)\n"
                  << "  Starting sort+count...\n";
        std::cout.flush();

        uint64_t hist[256] = {};
        uint64_t total_unique = 0;
        double   sum_sort = 0, sum_count = 0;

        std::vector<std::array<uint64_t, 256>> cu_hist(N_CU);
        for (int cu = 0; cu < N_CU; cu++) cu_hist[cu].fill(0);
        std::vector<uint64_t> cu_unique_v(N_CU, 0);
        std::vector<double>   cu_sort_t(N_CU, 0), cu_cnt_t(N_CU, 0);

        auto t_sort0 = std::chrono::high_resolution_clock::now();

        std::vector<std::thread> sort_threads;
        for (int cu = 0; cu < N_CU; cu++) {
            sort_threads.emplace_back([&, cu]() {
                /* Find this CU's max sub-bin size for temp buffer allocation */
                uint64_t max_sb = 0;
                for (int s = 0; s < N_SBINS; s++)
                    if (km_bins[cu][s].size() > max_sb)
                        max_sb = km_bins[cu][s].size();

                std::vector<uint64_t> tmp(max_sb);
                double sort_sum = 0, cnt_sum = 0;

                for (int s = 0; s < N_SBINS; s++) {
                    auto& sb  = km_bins[cu][s];
                    uint64_t n = sb.size();
                    if (n == 0) continue;

                    auto ts0 = std::chrono::high_resolution_clock::now();
                    lsd_radix_sort_42(sb.data(), tmp.data(), n);
                    auto ts1 = std::chrono::high_resolution_clock::now();
                    cu_unique_v[cu] += count_sorted(tmp.data(), n, cu_hist[cu].data());
                    auto ts2 = std::chrono::high_resolution_clock::now();
                    sort_sum += elapsed(ts0, ts1);
                    cnt_sum  += elapsed(ts1, ts2);
                }

                cu_sort_t[cu] = sort_sum;
                cu_cnt_t[cu]  = cnt_sum;
            });
        }
        for (auto& t : sort_threads) t.join();

        auto t_sort1 = std::chrono::high_resolution_clock::now();
        double sort_wall = elapsed(t_sort0, t_sort1);

        for (int cu = 0; cu < N_CU; cu++) {
            total_unique += cu_unique_v[cu];
            sum_sort     += cu_sort_t[cu];
            sum_count    += cu_cnt_t[cu];
            for (int i = 0; i < 256; i++) hist[i] += cu_hist[cu][i];
            std::cout << "  [CU " << cu << "] recs=" << [&]{
                uint64_t s=0; for(int i=0;i<N_SBINS;i++) s+=km_bins[cu][i].size(); return s;}()
                      << "  sort=" << cu_sort_t[cu] << "s"
                      << "  count=" << cu_cnt_t[cu] << "s"
                      << "  unique=" << cu_unique_v[cu] << "\n";
        }

        auto t_run_end = std::chrono::high_resolution_clock::now();
        double t_run = elapsed(t_run_start, t_run_end);

        std::cout << "\n  Batches        : " << batch_num << "\n"
                  << "  FASTQ total    : " << total_fastq_processed / 1048576.0 << " MB\n"
                  << "  Unique k-mers  : " << total_unique << "\n"
                  << "  Throughput     : " << total_fastq_processed / 1048576.0 / t_run << " MB/s\n"
                  << "  KA+DMA total   : " << sum_ka << " s\n"
                  << "  append total   : " << sum_append << " s\n"
                  << "  sort(CU sum)   : " << sum_sort << " s\n"
                  << "  count(CU sum)  : " << sum_count << " s\n"
                  << "  sort+count wall: " << sort_wall << " s\n"
                  << "  extraction e2e : " << t_extract << " s\n"
                  << "  e2e (all)      : " << t_run << " s\n";

        if (print_hist && run == n_runs - 1) {
            std::cout << "\n--- frequency histogram ---\n";
            for (int i = 1; i < 256; i++)
                if (hist[i]) std::cout << i << "\t" << hist[i] << "\n";
        }
        run_e2e_v.push_back(t_run);
    }

    if (n_runs > 1) {
        auto med = [](std::vector<double> v) {
            std::sort(v.begin(), v.end());
            int n = v.size();
            return n % 2 == 0 ? (v[n/2-1] + v[n/2]) / 2.0 : v[n/2];
        };
        std::cout << "\n[sort58] Median e2e over " << n_runs << " runs: "
                  << med(run_e2e_v) << " s\n";
    }
    return 0;
}
