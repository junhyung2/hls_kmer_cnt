/* exp73 — final host for the FPGA k-mer counter. 3 x Alveo U50, host-side counting.
 *
 * Runs on the SAME xclbin as exp50: every gain from exp50 -> exp73 is host-side,
 * ZERO resynthesis. The kernel (fastq_extractor_8port) streams FASTQ from gmem_fq,
 * extracts canonical k-mers, and writes h = hash64_fe(kmer) to one of 8 AXI ports
 * chosen by cu = h & 7. The host reads those back and counts them.
 *
 * Measured on ERR022075 (22.72M read pairs, k=21), 2026-07-10:
 *     e2e   10.25 s   (vs KMC 3.2.4 -t16: 14.21 s -> 1.39x)
 *     wall  13.2-14.1 s (vs KMC 14.1-14.3 s     -> 1.05x)
 *     unique 109,533,748 — IDENTICAL to KMC (exact, no k-mers lost)
 *
 * ── Shape ────────────────────────────────────────────────────────────────────
 *   3 devices x 3 CUs = 9 instances; each instance has 8 AXI k-mer ports (N_CU).
 *   Host hash table: N_CU x NBINS x BIN_SIZE = 8 x 4096 x 4096 = 134,217,728
 *   slots x 8 B = 1 GB (see the block above NBINS_BITS for the layout).
 *
 * ── The 3-stage software pipeline (the core of exp68) ────────────────────────
 *   Batch loop runs, for batch N:   kernel(N) || readback(N-1) || acc(N-2)
 *   km_bo/se_bo are double-buffered on [kb = N&1], so kernel(N) writes set kb
 *   while the host reads back set 1-kb (written by kernel(N-1)). flat_buf/
 *   scat_buf/flat_cnt are likewise 2-deep, so acc(N-2) reads a buffer disjoint
 *   from the one readback(N-1) is writing. Before exp68 the readback was
 *   serialized behind ext_run.wait(), which nothing required.
 *
 * ── Lineage (each step host-only, each measured) ─────────────────────────────
 *   exp50  batch-streamed baseline                        e2e 16.87 s
 *   exp65  KA/acc pipelining (acc(N) carried, not joined) e2e 12.29 s
 *   exp68  readback/kernel overlap -> 3-stage pipeline    e2e 11.65 s
 *   exp69  acc phase B parallelized over bin ranges (-as) e2e 10.95 s
 *   exp71a host splitmix64 rehash (see fe_bin/fe_slot)    e2e 10.27 s  + EXACT
 *   exp73  host buffers sized from batch_bytes, not from
 *          FE_MAX_PER_CU (was 6.4x over-allocated)        e2e 10.25 s, wall -2.2 s
 *
 * ── Where the time goes (exp67/exp75 profiling; NOT a guess) ─────────────────
 *   kernel exec  ~5.2 s  — bound by the gmem_fq AXI read (1.23 B/cycle), NOT by
 *                          pipeline width. exp66 (dual-lane kernel) therefore
 *                          gave zero e2e gain despite being correct on silicon.
 *   readback     ~6.5 s  (29.3 GB D2H, 3.66G k-mers)
 *   upload       ~6.4 s  — fully hidden inside pf_th, not on the critical path.
 *   acc          dominates e2e. Inside acc: scatter 45.5 core-s vs insert 19.0
 *                core-s — the SCATTER is the bottleneck, not the insert.
 *                (exp75 measured this; earlier comments here claimed the opposite.)
 *
 * Usage: h73_slim <xclbin> <fastq1> <fastq2> <k>
 *                 [-b <mb>] [-n <runs>] [-as <n>] [-hist]
 *   -b  <mb>    max FASTQ per instance per batch in MiB (default: whole file, 1 batch)
 *   -n  <runs>  repeat the full multi-batch sweep N times (reports median e2e)
 *   -as <n>     acc phase-B sub-threads per cu; N_CU*n threads during insert
 *   -hist       print the frequency histogram after the last run
 *
 * Champion config:  h73_slim <xclbin> <fq1> <fq2> 21 -b 40 -n 3 -as 4
 * Build:            -DUSE_REHASH=1 -DBIN_BITS_CFG=12
 *
 * Correctness check: unique MUST be 109,533,748 on ERR022075 (KMC-exact) and
 * "dropped inserts" MUST be 0, for any -b / -as. Do NOT use the pre-exp71
 * numbers (e.g. exp49/exp50's 105,412,284) as the reference — those were the
 * result of the silent bin-overflow drops that exp71 fixed.
 */

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <thread>
#include <atomic>
#include <stdexcept>
#include <cstdio>
#include <cassert>

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
static constexpr int N_CU           = FE_N_CU;  /* 8 */

/* ── Persistent accumulation hash table (2-level, cache-friendly) ──────────
 * One table per CU. Each table is NBINS bins of BIN_SIZE open-addressed slots:
 *
 *   NBINS_BITS = 12 -> NBINS    = 4096 bins
 *   BIN_BITS   = 12 -> BIN_SIZE = 4096 slots       (BIN_BITS_CFG, build-time)
 *
 *   per CU : 4096 x 4096 = 16,777,216 slots x 8 B = 128 MB
 *   total  : 8 CU x 128 MB                        =   1 GB host RAM
 *   per bin: 4096 x 8 B = 32 KB -> L1-resident during insert  <- the whole point
 *
 * The 32 KB/bin figure is why acc scatters into bins first: the random-probe
 * insert then walks a 32 KB table that stays in L1 (exp75 confirmed the insert
 * is already at L1 speed, and that adding an L3-sized hot cache in front of it
 * made things WORSE).
 *
 * A slot is one uint64: (count << KM_BITS) | kmer42.  EMPTY = 0.
 *
 * Bin/slot are derived from h by fe_bin()/fe_slot() below — see the exp71 note
 * there. They are NOT plain bit-slices of h; that was the pre-exp71 scheme and
 * it silently dropped k-mers.
 *
 * Bins are NEVER cleared inside a run -> counts accumulate across all batches,
 * so memory stays O(unique k-mers), not O(total FASTQ).
 *
 * Occupancy on ERR022075 (109.5M unique over 134.2M slots): ~81.6% load with
 * the rehash, and 0 bins full. Without it: 78.5% load but 18% of bins 100% full
 * -> 4.35M dropped inserts. Load alone does not tell you it is safe; see exp71.
 */
static constexpr int      NBINS_BITS  = 12;
static constexpr int      NBINS       = 1 << NBINS_BITS;
#ifndef BIN_BITS_CFG
#define BIN_BITS_CFG 12
#endif
static constexpr int      BIN_BITS    = BIN_BITS_CFG;
static constexpr int      BIN_SIZE    = 1 << BIN_BITS;   /* 4096 */
static constexpr uint64_t BIN_MASK    = BIN_SIZE - 1;
static constexpr int      KM_BITS     = 42;
static constexpr uint64_t KEY_MASK    = (1ULL << KM_BITS) - 1;
static constexpr uint64_t ACC_CNT_MAX = (1ULL << (64 - KM_BITS)) - 1;
static constexpr uint64_t EMPTY       = 0ULL;
static std::atomic<uint64_t> g_drops{0};   /* inserts discarded: bin full */

/* exp71: bin/slot derivation.
 * The FPGA sends h = hash64_fe(canonical_kmer), which is a BIJECTION on 42 bits
 * (GF(2)-linear, both shifts nilpotent) -> it carries full k-mer identity but has
 * NO avalanche. Deriving bin/slot directly from h's bit fields therefore inherits
 * genomic structure and loads bins very unevenly (measured: 18% of bins 100% full
 * at only 78.5% global load -> 4.35M inserts silently dropped).
 *
 * cu = h&7 is fixed by which AXI port the k-mer arrived on (measured balanced,
 * +-1.6%), but bin/slot are chosen HERE, host-side. So we can re-mix h with a
 * strong non-linear finalizer (splitmix64) with NO kernel change / no resynthesis.
 * The stored key remains the full 42-bit h, so identity/uniqueness is preserved. */
#ifndef USE_REHASH
#define USE_REHASH 1
#endif

static inline uint64_t fe_mix64(uint64_t x)
{
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

#if USE_REHASH
static inline uint64_t fe_bin (uint64_t h) { return fe_mix64(h) & (NBINS - 1); }
static inline uint32_t fe_slot(uint64_t h) { return (uint32_t)((fe_mix64(h) >> NBINS_BITS) & BIN_MASK); }
#else
static inline uint64_t fe_bin (uint64_t h) { return (h >> FE_N_CU_BITS) & (NBINS - 1); }
static inline uint32_t fe_slot(uint64_t h) { return (uint32_t)((h >> (FE_N_CU_BITS + NBINS_BITS)) & BIN_MASK); }
#endif

/* ── mmap wrapper ─────────────────────────────────────────────────────────*/
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

/* Scan at most ~2 KB to find '\n@' at or after target — O(1) per call */
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

/* ── acc: fold a batch of k-mers into the persistent per-CU table ────────────
 * Two phases, split in exp69 so the second can be parallelized:
 *
 *   phase A  acc_scatter_phase() — count keys per bin, prefix-sum, then scatter
 *            them into contiguous per-bin groups.  1 thread per cu.
 *   phase B  acc_insert_bins()   — per bin, open-address insert/increment into
 *            the persistent HT.  Split across `acc_split` disjoint bin ranges.
 *
 * Why phase B is safe to parallelize: bin b owns exactly
 * pers_flat[b*BIN_SIZE .. (b+1)*BIN_SIZE) and scatter has already grouped keys
 * by bin, so disjoint bin ranges touch disjoint HT memory and disjoint scatter
 * memory -> race-free. The resulting HT state and new_unique total are identical
 * to a serial run (verified: unique is stable across -as 1/2/4).
 *
 * COST (exp75, core-seconds summed over threads):  A = 45.5,  B = 19.0.
 * Phase A (scatter) is the real bottleneck — 2.4x phase B. Earlier revisions of
 * this file asserted the insert was dominant; that was wrong. Attacking A
 * directly still failed: exp76 (2-pass radix) made it worse because the 4096
 * open streams x 64 B already fit in L2, and exp77 (parallel scatter) cut acc
 * 12.5% but did not move e2e, because acc is no longer the sole critical path.
 *
 * phase A keeps its stack-local bin_count/bin_start/pos arrays (exp64: turning
 * those into pointer params costs ~15% via lost aliasing info); the results are
 * copied out once at the end (2 x 32 KB, negligible).
 *
 * pers_flat: flat NBINS*BIN_SIZE entries, bin b starts at b*BIN_SIZE.
 *            Persistent across batches — never cleared inside a run.
 * scatter:   temporary work buffer, size >= n.
 * acc_insert_bins() returns the number of NEWLY inserted unique k-mers.       */
static void acc_scatter_phase(const uint64_t* src, uint64_t* scatter, uint64_t n,
                              uint64_t* out_bin_count, uint64_t* out_bin_start)
{
    uint64_t bin_count[NBINS] = {};
    for (uint64_t j = 0; j < n; j++)
        bin_count[fe_bin(src[j])]++;

    uint64_t bin_start[NBINS + 1];
    bin_start[0] = 0;
    for (int i = 0; i < NBINS; i++)
        bin_start[i + 1] = bin_start[i] + bin_count[i];

    uint64_t pos[NBINS];
    std::copy(bin_start, bin_start + NBINS, pos);
    for (uint64_t j = 0; j < n; j++) {
        uint64_t v = src[j];
        scatter[pos[fe_bin(v)]++] = v;
    }

    std::memcpy(out_bin_count, bin_count, NBINS * sizeof(uint64_t));
    std::memcpy(out_bin_start, bin_start, NBINS * sizeof(uint64_t));
}

/* phase B: insert bins [b_lo, b_hi). Called once per acc_split sub-range. */
static uint64_t acc_insert_bins(const uint64_t* scatter,
                                const uint64_t* bin_count, const uint64_t* bin_start,
                                uint64_t* pers_flat, int b_lo, int b_hi)
{
    uint64_t new_unique = 0;
    for (int b = b_lo; b < b_hi; b++) {
        uint64_t b_n = bin_count[b];
        if (b_n == 0) continue;

        uint64_t* ht = pers_flat + (uint64_t)b * BIN_SIZE;
        const uint64_t* bdata = scatter + bin_start[b];
        for (uint64_t j = 0; j < b_n; j++) {
            uint64_t kmer = bdata[j] & KEY_MASK;
            uint32_t h = fe_slot(bdata[j]);
            uint32_t start_h = h;
            while (ht[h] != EMPTY && (ht[h] & KEY_MASK) != kmer) {
                h = (h + 1) & BIN_MASK;
                if (h == start_h) { h = BIN_SIZE; break; }  /* bin full — drop */
            }
            if (h == (uint32_t)BIN_SIZE) { g_drops.fetch_add(1, std::memory_order_relaxed); continue; }
            if (ht[h] == EMPTY) {
                ht[h] = (1ULL << KM_BITS) | kmer;
                new_unique++;
            } else {
                uint64_t cnt = ht[h] >> KM_BITS;
                if (cnt < ACC_CNT_MAX)
                    ht[h] += (1ULL << KM_BITS);
            }
        }
    }
    return new_unique;
}

/* Scan all bins in pers_flat (NBINS*BIN_SIZE), count unique entries, fill histogram. */
static uint64_t finalize_acc(const uint64_t* pers_flat, uint64_t* histogram)
{
    uint64_t unique = 0;
    for (int b = 0; b < NBINS; b++) {
        const uint64_t* ht = pers_flat + (uint64_t)b * BIN_SIZE;
        for (int s = 0; s < BIN_SIZE; s++) {
            if (ht[s] == EMPTY) continue;
            unique++;
            uint64_t cnt = ht[s] >> KM_BITS;
            histogram[cnt < 255 ? cnt : 255]++;
        }
    }
    return unique;
}

int main(int argc, char* argv[])
{
    /* exp72: per-phase timing of everything OUTSIDE the batch-loop e2e */
    auto T_prog = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string,double>> g_phase;
    auto PH = [&](const char* n, std::chrono::high_resolution_clock::time_point a) {
        g_phase.emplace_back(std::string(n), elapsed(a, std::chrono::high_resolution_clock::now()));
    };
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <xclbin> <fastq1> <fastq2> <k>"
                     " [-b <mb_per_inst>] [-n <runs>] [-as <n>] [-hist]\n";
        return 1;
    }
    std::string xclbin_path = argv[1];
    std::string fq1_path    = argv[2];
    std::string fq2_path    = argv[3];
    int         k           = std::stoi(argv[4]);
    uint64_t    batch_mb    = 0;   /* 0 = auto (single batch) */
    int         n_runs      = 1;
    bool        print_hist  = false;
    int         acc_split   = 2;   /* H2: bin-range sub-threads per cu (8*acc_split total) */
    for (int i = 5; i < argc; i++) {
        if (std::string(argv[i]) == "-b"    && i+1 < argc) batch_mb   = std::stoull(argv[++i]);
        if (std::string(argv[i]) == "-n"    && i+1 < argc) n_runs     = std::stoi(argv[++i]);
        if (std::string(argv[i]) == "-hist") print_hist = true;
        if (std::string(argv[i]) == "-as"   && i+1 < argc) acc_split  = std::stoi(argv[++i]);
    }
    if (k != FE_K) { std::cerr << "Only k=" << FE_K << " supported\n"; return 1; }
    if (acc_split < 1 || NBINS % acc_split != 0) {
        std::cerr << "-as must divide NBINS(" << NBINS << ")\n"; return 1;
    }
    std::cout << "  REHASH=" << USE_REHASH << "  BIN_BITS=" << BIN_BITS
              << "  capacity=" << ((uint64_t)N_CU*NBINS*BIN_SIZE) << " slots\n";
    std::cout << "  acc_split      : " << acc_split
              << "  (" << N_CU * acc_split << " acc threads during insert)\n";

    /* ── 1. mmap both files ─────────────────────────────────────────────── */
    MmapFile fq1(fq1_path), fq2(fq2_path);
    uint64_t fq_total = fq1.size + fq2.size;
    PH("argparse + FASTQ open/size probe", T_prog);
    std::cout << "[stream] Total FASTQ: " << fq_total / 1048576.0 << " MB\n";

    /* batch_bytes = per-instance max (auto → entire file in one batch) */
    uint64_t batch_bytes = (batch_mb > 0)
        ? batch_mb * 1048576ULL
        : (fq_total / N_INST_TOTAL + 4096);  /* single batch */

    /* Estimate max batches to report */
    int est_batches = (int)((fq_total + (uint64_t)N_INST_TOTAL * batch_bytes - 1)
                            / ((uint64_t)N_INST_TOTAL * batch_bytes));
    std::cout << "  batch_bytes/inst: " << batch_bytes / 1048576.0 << " MB"
              << "  estimated batches: " << est_batches << "\n";

    /* ── 2. Open devices ─────────────────────────────────────────────────── */
    xrt::device device[N_DEV];
    xrt::uuid   uuid[N_DEV];
    {
        auto a = std::chrono::high_resolution_clock::now();
        for (int d = 0; d < N_DEV; d++) device[d] = xrt::device(d);
        PH("xrt::device open (3 cards)", a);
        a = std::chrono::high_resolution_clock::now();
        for (int d = 0; d < N_DEV; d++) uuid[d] = device[d].load_xclbin(xclbin_path);
        PH("load_xclbin (3 cards)", a);
    }

    const char* krnl_names[N_INST_PER_DEV] = {
        "fastq_extractor_8port:{fastq_extractor_8port_1}",
        "fastq_extractor_8port:{fastq_extractor_8port_2}",
        "fastq_extractor_8port:{fastq_extractor_8port_3}",
    };

    /* ── 3. Pre-allocate BOs at fixed batch_bytes ────────────────────────── */
    const uint64_t cu_bytes = (uint64_t)FE_MAX_PER_CU * sizeof(uint64_t);
    const uint64_t se_bytes = (uint64_t)N_CU * sizeof(uint64_t);

    xrt::kernel ext_krnl[N_DEV][N_INST_PER_DEV];
    xrt::bo     fq_bo   [N_DEV][N_INST_PER_DEV][2]; /* ping-pong for double-buffer */
    /* exp68 (H1): km_bo/se_bo double-buffered so kernel(N+1) can write set [kb]
     * while the host reads back set [1-kb] from batch N. Each km bank then holds
     * 2 x 64MB = 128MB, well inside the 256MB HBM bank size.                    */
    xrt::bo     km_bo   [N_DEV][N_INST_PER_DEV][N_CU][2];
    xrt::bo     se_bo   [N_DEV][N_INST_PER_DEV][2];

    /* BO for fq gets a small extra margin so that FASTQ record boundary overshoot
     * (find_next_header_at can advance up to ~400 bytes past the target) never
     * causes the per-instance sub_len to exceed the allocated BO size.        */
    const uint64_t fq_bo_bytes = (batch_bytes + 8192 + 7) & ~7ULL;

    auto t_k0 = std::chrono::high_resolution_clock::now();
    for (int d = 0; d < N_DEV; d++)
        for (int ci = 0; ci < N_INST_PER_DEV; ci++)
            ext_krnl[d][ci] = xrt::kernel(device[d], uuid[d], krnl_names[ci]);
    PH("xrt::kernel objects (9 CUs)", t_k0);

    auto t_bo0 = std::chrono::high_resolution_clock::now();
    for (int d = 0; d < N_DEV; d++) {
        for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
            for (int pp = 0; pp < 2; pp++)
                fq_bo[d][ci][pp] = xrt::bo(device[d], fq_bo_bytes,
                                            ext_krnl[d][ci].group_id(0));
            for (int kb = 0; kb < 2; kb++)
                se_bo[d][ci][kb] = xrt::bo(device[d], se_bytes,
                                            ext_krnl[d][ci].group_id(10));
            for (int cu = 0; cu < N_CU; cu++)
                for (int kb = 0; kb < 2; kb++)
                    km_bo[d][ci][cu][kb] = xrt::bo(device[d], cu_bytes,
                                                    ext_krnl[d][ci].group_id(2 + cu));
        }
    }
    PH("BO alloc (fq x2, se x2, km 8x2 -> HBM)", t_bo0);

    /* ── 4. Host-side buffers ────────────────────────────────────────────── */
    /* Double-buffered (buf 0/1): inst_workers(N+1) write into "cur" while
     * acc_workers(N) still read "prev" — same ping-pong idea as fq_bo.      */
    /* exp73: HOST buffers were strided by FE_MAX_PER_CU (device km_bo capacity,
     * 8.39M entries = 64MB per inst per cu) -> 2 buf x 8 cu x 2 arrays x 9 inst
     * x 64MB = 19.3 GB, all value-initialized to 0 at startup (measured 2.36 s).
     * But a batch only ever produces ~0.25 k-mers/FASTQ-byte, spread over 8 cu:
     *   real per (inst,cu) per batch ~= 1.30M entries  -> 6.4x over-allocated.
     * Size the host stride from batch_bytes instead. km_bo (device) is left at
     * FE_MAX_PER_CU: the kernel writes it with no bound check, so shrinking that
     * could silently corrupt HBM. The host side is guarded at runtime below.    */
    const uint64_t HOST_STRIDE = std::min<uint64_t>(
        (uint64_t)FE_MAX_PER_CU,
        std::max<uint64_t>(1u << 20, fq_bo_bytes / 16));   /* ~2x real need */
    const uint64_t max_per_cu_total = HOST_STRIDE * N_INST_TOTAL;
    std::atomic<uint64_t> g_max_n{0};   /* largest observed n per (inst,cu) */
    std::cout << "  HOST_STRIDE    : " << HOST_STRIDE << " entries/inst/cu  ("
              << (max_per_cu_total * 8.0 / 1048576.0) << " MB per (buf,cu), total "
              << (2.0 * N_CU * 2 * max_per_cu_total * 8.0 / (1<<30)) << " GB)\n";
    std::vector<uint64_t> flat_buf[2][N_CU];
    std::vector<uint64_t> scat_buf[2][N_CU];  /* scatter work buffer for 2-level acc */
    {
      auto a = std::chrono::high_resolution_clock::now();
      for (int buf = 0; buf < 2; buf++)
        for (int cu = 0; cu < N_CU; cu++) {
            flat_buf[buf][cu].resize(max_per_cu_total, 0);
            scat_buf[buf][cu].resize(max_per_cu_total, 0);
        }
      PH("host flat_buf+scat_buf alloc+zero", a);
    }

    uint64_t flat_cnt[2][N_INST_TOTAL][N_CU] = {};

    /* exp69 (H2): per-cu bin metadata shared between acc phase A (scatter) and
     * phase B (insert). Only ONE acc generation is ever live — acc(N-2) is joined
     * before acc(N-1) is spawned — so a single set per cu is safe.            */
    std::vector<uint64_t> bin_cnt_buf[N_CU];
    std::vector<uint64_t> bin_str_buf[N_CU];
    for (int cu = 0; cu < N_CU; cu++) {
        bin_cnt_buf[cu].resize(NBINS, 0);
        bin_str_buf[cu].resize(NBINS, 0);
    }

    /* ── 5. Persistent 2-level accumulation tables (per-CU, per-run) ────── */
    /* pers_flat[cu]: NBINS*BIN_SIZE = 4096*4096 = 16.78M entries = 128 MB/CU */
    /* Total: 8 CU x 128 MB = 1 GB. Layout/rationale: see the NBINS block.    */
    const uint64_t PERS_CU_SIZE = (uint64_t)NBINS * BIN_SIZE;
    std::vector<uint64_t> pers_flat[N_CU];
    {
      auto a = std::chrono::high_resolution_clock::now();
      for (int cu = 0; cu < N_CU; cu++)
        pers_flat[cu].assign(PERS_CU_SIZE, EMPTY);
      PH("pers_flat (HT) alloc+init", a);
    }

    /* ── 6. Multi-run loop ───────────────────────────────────────────────── */
    std::vector<double> run_e2e_v;

    for (int run = 0; run < n_runs; run++) {
        std::cout << "\n[stream] ====== Run " << (run+1) << " / " << n_runs << " ======\n";

        /* Reset persistent bins for this run */
        {
          auto a = std::chrono::high_resolution_clock::now();
          for (int cu = 0; cu < N_CU; cu++)
            std::fill(pers_flat[cu].begin(), pers_flat[cu].end(), EMPTY);
          PH("HT reset std::fill (per run)", a);
        }

        uint64_t total_unique = 0;
        uint64_t total_fastq_processed = 0;
        int batch_num = 0;

        double sum_ka = 0, sum_acc = 0;

        /* ── Batch helpers: double-buffer (fill+PCIe of N+1 overlaps KA of N) ─ */
        struct BatchInfo {
            uint64_t start = 0, end = 0, size = 0;
            uint64_t sub_split[N_INST_TOTAL + 1]{};
            uint64_t sub_len[N_INST_TOTAL]{};
            uint64_t sub_pad[N_INST_TOTAL]{};
        };

        auto compute_batch_info = [&](uint64_t bstart) -> BatchInfo {
            BatchInfo bi;
            bi.start = bstart;
            uint64_t desired = bstart + (uint64_t)N_INST_TOTAL * batch_bytes;
            bi.end = (desired >= fq_total) ? fq_total
                     : find_next_header_at(fq1, fq2, fq_total, desired);
            if (bi.end == fq_total && desired < fq_total) bi.end = fq_total;
            bi.size = bi.end - bi.start;
            bi.sub_split[0] = bstart;
            for (int i = 1; i < N_INST_TOTAL; i++) {
                uint64_t t = bstart + (uint64_t)i * bi.size / N_INST_TOTAL;
                bi.sub_split[i] = find_next_header_at(fq1, fq2, fq_total, t);
            }
            bi.sub_split[N_INST_TOTAL] = bi.end;
            for (int i = 0; i < N_INST_TOTAL; i++) {
                bi.sub_len[i] = bi.sub_split[i+1] - bi.sub_split[i];
                bi.sub_pad[i] = (bi.sub_len[i] + 7) & ~7ULL;
            }
            return bi;
        };

        /* Fill fq_bo[pp] from mmap + 9-parallel PCIe sync to device */
        auto fill_and_sync_pp = [&](const BatchInfo& bi, int pp) {
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
                    int inst = d * N_INST_PER_DEV + ci;
                    if (bi.sub_pad[inst] > fq_bo_bytes) {
                        std::cerr << "[ERROR] inst " << inst << " sub_pad="
                                  << bi.sub_pad[inst] << " > fq_bo_bytes="
                                  << fq_bo_bytes << "\n";
                        std::exit(1);
                    }
                    uint8_t* dst = fq_bo[d][ci][pp].map<uint8_t*>();
                    cat_memcpy(dst, fq1, fq2, bi.sub_split[inst], bi.sub_len[inst]);
                    if (bi.sub_pad[inst] > bi.sub_len[inst])
                        std::memset(dst + bi.sub_len[inst], 0,
                                    bi.sub_pad[inst] - bi.sub_len[inst]);
                }
            std::vector<std::thread> up_th;
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
                    int inst = d * N_INST_PER_DEV + ci;
                    up_th.emplace_back([&, inst, pp]() {
                        int d2 = inst / N_INST_PER_DEV, ci2 = inst % N_INST_PER_DEV;
                        fq_bo[d2][ci2][pp].sync(XCL_BO_SYNC_BO_TO_DEVICE,
                                                 bi.sub_pad[inst], 0);
                    });
                }
            for (auto& t : up_th) t.join();
        };

        auto t_run_start = std::chrono::high_resolution_clock::now();

        /* Pre-fill batch 0 into pp=0 (serial, included in t_run) */
        BatchInfo cur_bi = compute_batch_info(0);
        fill_and_sync_pp(cur_bi, 0);
        int pp_cur  = 0;

        /* exp68 (H1): 3-stage software pipeline
         *     kernel(N)  ||  readback(N-1)  ||  acc(N-2)
         * km_bo/se_bo are double-buffered on [kb], so kernel(N) writes set kb=N&1
         * while the host reads back set 1-kb (written by kernel(N-1)).
         * flat_buf/flat_cnt/scat_buf stay 2-deep:
         *     readback(N-1) writes flat_buf[1-kb]   ((N-1)%2)
         *     acc(N-2)      reads  flat_buf[kb]     ((N-2)%2)
         * which are always disjoint.                                            */
        std::vector<std::thread> acc_workers_prev;
        std::atomic<uint64_t>    new_unique_prev{0};
        bool                     have_prev = false;
        double sum_acc_wait = 0, sum_rb_wait = 0, sum_k_wait = 0, sum_upload = 0;
        uint64_t total_kmers = 0;

        /* readback(batch written into km_bo set kb_src) -> flat_buf[buf_dst] */
        auto spawn_readback = [&](int kb_src, int buf_dst, std::vector<std::thread>& ws) {
            std::fill(&flat_cnt[buf_dst][0][0],
                      &flat_cnt[buf_dst][0][0] + N_INST_TOTAL * N_CU, 0ULL);
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
                    ws.emplace_back([&, d, ci, kb_src, buf_dst]() {
                        int inst = d * N_INST_PER_DEV + ci;
                        se_bo[d][ci][kb_src].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
                        const uint64_t* se = se_bo[d][ci][kb_src].map<const uint64_t*>();
                        for (int cu = 0; cu < N_CU; cu++)
                            flat_cnt[buf_dst][inst][cu] = se[cu];

                        for (int cu = 0; cu < N_CU; cu++) {
                            uint64_t n = flat_cnt[buf_dst][inst][cu];
                            if (n > 0) {
                                km_bo[d][ci][cu][kb_src].sync(XCL_BO_SYNC_BO_FROM_DEVICE,
                                                               n * sizeof(uint64_t), 0);
                                const uint64_t* src =
                                    km_bo[d][ci][cu][kb_src].map<const uint64_t*>();
                                if (n > HOST_STRIDE) {
                                    fprintf(stderr,
                                        "FATAL exp73: inst %d cu %d n=%llu > HOST_STRIDE=%llu"
                                        " (raise stride)\n", inst, cu,
                                        (unsigned long long)n,
                                        (unsigned long long)HOST_STRIDE);
                                    std::abort();
                                }
                                {   uint64_t prev = g_max_n.load(std::memory_order_relaxed);
                                    while (n > prev && !g_max_n.compare_exchange_weak(prev, n)) {}
                                }
                                uint64_t dst_off = (uint64_t)inst * HOST_STRIDE;
                                std::memcpy(flat_buf[buf_dst][cu].data() + dst_off,
                                            src, n * sizeof(uint64_t));
                            }
                        }
                    });
                }
        };

        /* acc over flat_buf[buf_src] — spawned here, joined in a later iteration.
         * exp69: phase A (compact+scatter) stays 1 thread/cu; phase B (the
         * random-probe insert) is split across `acc_split` disjoint bin ranges,
         * giving N_CU*acc_split threads during phase B.
         * NB: phase A is the more expensive half (45.5 vs 19.0 core-s, exp75) and
         * is still single-threaded per cu. exp77 parallelized it — acc dropped
         * 12.5% but e2e did not move, so it was not merged. This is the obvious
         * place to look if acc ever becomes the sole critical path again.      */
        auto spawn_acc = [&](int buf_src) {
            new_unique_prev.store(0, std::memory_order_relaxed);
            acc_workers_prev.clear();
            for (int cu = 0; cu < N_CU; cu++) {
                acc_workers_prev.emplace_back([&, cu, buf_src]() {
                    uint64_t write_off = 0;
                    for (int inst = 0; inst < N_INST_TOTAL; inst++) {
                        uint64_t n = flat_cnt[buf_src][inst][cu];
                        if (n == 0) continue;
                        uint64_t src_off = (uint64_t)inst * HOST_STRIDE;
                        if (write_off != src_off)
                            std::memmove(flat_buf[buf_src][cu].data() + write_off,
                                         flat_buf[buf_src][cu].data() + src_off,
                                         n * sizeof(uint64_t));
                        write_off += n;
                    }
                    if (write_off == 0) return;

                    uint64_t*       bc   = bin_cnt_buf[cu].data();
                    uint64_t*       bs   = bin_str_buf[cu].data();
                    uint64_t*       scat = scat_buf[buf_src][cu].data();
                    uint64_t*       ht   = pers_flat[cu].data();

                    /* phase A */
                    acc_scatter_phase(flat_buf[buf_src][cu].data(), scat, write_off, bc, bs);

                    /* phase B: disjoint bin ranges -> disjoint HT + scatter memory */
                    const int per = NBINS / acc_split;
                    std::vector<uint64_t>    nu(acc_split, 0);
                    std::vector<std::thread> helpers;
                    for (int s = 1; s < acc_split; s++)
                        helpers.emplace_back([&, s]() {
                            nu[s] = acc_insert_bins(scat, bc, bs, ht, s * per, (s + 1) * per);
                        });
                    nu[0] = acc_insert_bins(scat, bc, bs, ht, 0, per);
                    for (auto& h : helpers) h.join();

                    uint64_t tot = 0;
                    for (int s = 0; s < acc_split; s++) tot += nu[s];
                    new_unique_prev.fetch_add(tot, std::memory_order_relaxed);
                });
            }
            have_prev = true;
        };

        /* ── Batch loop (3-stage: kernel(N) || readback(N-1) || acc(N-2)) ─── */
        int  N        = 0;
        int  last_kb  = 0;
        while (cur_bi.start < fq_total) {
            batch_num++;
            uint64_t batch_size = cur_bi.size;
            total_fastq_processed += batch_size;

            const int kb = N & 1;          /* kernel(N) writes km_bo set kb        */
            last_kb = kb;

            auto t_b0 = std::chrono::high_resolution_clock::now();

            /* Prefetch FASTQ for batch N+1 (runs under everything else) */
            int pp_nxt = 1 - pp_cur;
            BatchInfo nxt_bi;
            bool has_next = (cur_bi.end < fq_total);
            std::thread pf_th;
            double pf_time = 0;
            if (has_next) {
                pf_th = std::thread([&, pp_nxt]() {
                    auto p0 = std::chrono::high_resolution_clock::now();
                    nxt_bi = compute_batch_info(cur_bi.end);
                    fill_and_sync_pp(nxt_bi, pp_nxt);
                    pf_time = elapsed(p0, std::chrono::high_resolution_clock::now());
                });
            }

            /* Stage A: launch kernel(N) into km_bo/se_bo set kb (async) */
            xrt::run ext_run[N_DEV][N_INST_PER_DEV];
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++) {
                    int inst = d * N_INST_PER_DEV + ci;
                    ext_run[d][ci] = ext_krnl[d][ci](
                        fq_bo[d][ci][pp_cur], cur_bi.sub_len[inst],
                        km_bo[d][ci][0][kb], km_bo[d][ci][1][kb],
                        km_bo[d][ci][2][kb], km_bo[d][ci][3][kb],
                        km_bo[d][ci][4][kb], km_bo[d][ci][5][kb],
                        km_bo[d][ci][6][kb], km_bo[d][ci][7][kb],
                        se_bo[d][ci][kb]);
                }

            /* Stage B: readback(N-1) from set 1-kb -> flat_buf[1-kb], concurrent
             * with kernel(N). Safe: kernel(N-1) was joined at the end of iter N-1. */
            std::vector<std::thread> rb_workers;
            if (N >= 1) spawn_readback(1 - kb, 1 - kb, rb_workers);

            /* Stage C: join acc(N-2) (reads flat_buf[kb], disjoint from readback) */
            if (have_prev) {
                auto tw0 = std::chrono::high_resolution_clock::now();
                for (auto& t : acc_workers_prev) t.join();
                sum_acc_wait += elapsed(tw0, std::chrono::high_resolution_clock::now());
                total_unique += new_unique_prev.load();
                have_prev = false;
            }

            /* Join readback(N-1), then hand its buffer to acc(N-1) */
            double rb_w = 0;
            if (N >= 1) {
                auto tr0 = std::chrono::high_resolution_clock::now();
                for (auto& t : rb_workers) t.join();
                rb_w = elapsed(tr0, std::chrono::high_resolution_clock::now());
                sum_rb_wait += rb_w;
                for (int i = 0; i < N_INST_TOTAL; i++)
                    for (int cu = 0; cu < N_CU; cu++)
                        total_kmers += flat_cnt[1 - kb][i][cu];
                spawn_acc(1 - kb);
            }

            /* Join kernel(N): its km_bo set kb must be complete before iter N+1
             * reads it back.                                                    */
            auto tk0 = std::chrono::high_resolution_clock::now();
            for (int d = 0; d < N_DEV; d++)
                for (int ci = 0; ci < N_INST_PER_DEV; ci++)
                    ext_run[d][ci].wait();
            double k_w = elapsed(tk0, std::chrono::high_resolution_clock::now());
            sum_k_wait += k_w;

            if (pf_th.joinable()) pf_th.join();
            sum_upload += pf_time;
            sum_ka += elapsed(t_b0, std::chrono::high_resolution_clock::now());

            std::cout << "  [batch " << batch_num << "]"
                      << "  FASTQ: " << batch_size / 1048576.0 << " MB"
                      << "  iter: " << elapsed(t_b0, std::chrono::high_resolution_clock::now()) << "s"
                      << "  | k_wait: " << k_w << "s"
                      << "  rb_wait: " << rb_w << "s"
                      << "  acc_wait(cum): " << sum_acc_wait << "s"
                      << "\n";

            if (!has_next) break;
            pp_cur = pp_nxt;
            cur_bi = nxt_bi;
            N++;
        }  /* end batch loop */

        /* ── Drain: kernel(M) already joined. Still owe readback(M) + acc(M-1,M) ── */
        {
            const int kb_M = last_kb;
            auto td0 = std::chrono::high_resolution_clock::now();

            if (have_prev) {                 /* acc(M-1), reads flat_buf[1-kb_M] */
                for (auto& t : acc_workers_prev) t.join();
                total_unique += new_unique_prev.load();
                have_prev = false;
            }

            std::vector<std::thread> rb_workers;   /* readback(M) -> flat_buf[kb_M] */
            spawn_readback(kb_M, kb_M, rb_workers);
            for (auto& t : rb_workers) t.join();
            for (int i = 0; i < N_INST_TOTAL; i++)
                for (int cu = 0; cu < N_CU; cu++)
                    total_kmers += flat_cnt[kb_M][i][cu];

            spawn_acc(kb_M);                        /* acc(M) */
            for (auto& t : acc_workers_prev) t.join();
            total_unique += new_unique_prev.load();
            have_prev = false;

            sum_acc += elapsed(td0, std::chrono::high_resolution_clock::now());
        }
        sum_acc += sum_acc_wait;

        auto t_run_end = std::chrono::high_resolution_clock::now();
        double t_run = elapsed(t_run_start, t_run_end);

        /* Verify total unique from accumulated table */
        uint64_t hist[256] = {};
        uint64_t verified_unique = 0;
        {
          auto a = std::chrono::high_resolution_clock::now();
          for (int cu = 0; cu < N_CU; cu++)
            verified_unique += finalize_acc(pers_flat[cu].data(), hist);
          PH("verify scan finalize_acc (per run)", a);
        }

        std::cout << "\n  Batches        : " << batch_num << "\n";
        std::cout << "  FASTQ total    : " << total_fastq_processed / 1048576.0 << " MB\n";
        std::cout << "  Unique (acc)   : " << total_unique << "  (verified: " << verified_unique << ")\n";
        std::cout << "  Throughput     : " << total_fastq_processed / 1048576.0 / t_run << " MB/s\n";
        std::cout << "  KA+DMA total   : " << sum_ka  << " s  (fill+PCIe overlapped)\n";
        std::cout << "  acc total      : " << sum_acc << " s\n";
        std::cout << "  e2e (all batch): " << t_run   << " s\n";

        /* ── exp68 (H1): blocking time per pipeline stage (sums ~= e2e) ── */
        double rb_gb = (double)total_kmers * sizeof(uint64_t) / 1e9;
        std::cout << "\n  --- 3-stage pipeline blocking time (kernel||readback||acc) ---\n";
        std::cout << "  kernel wait    : " << sum_k_wait   << " s\n";
        std::cout << "  readback wait  : " << sum_rb_wait  << " s   ("
                  << rb_gb << " GB D2H)\n";
        std::cout << "  acc wait       : " << sum_acc_wait << " s\n";
        std::cout << "  upload (hidden): " << sum_upload   << " s\n";
        std::cout << "  k-mers D2H     : " << total_kmers  << "\n";
        {
            uint64_t full_bins = 0, occupied = 0;
            for (int cu = 0; cu < N_CU; cu++)
                for (int b = 0; b < NBINS; b++) {
                    const uint64_t* ht = pers_flat[cu].data() + (uint64_t)b * BIN_SIZE;
                    uint64_t occ = 0;
                    for (int s2 = 0; s2 < BIN_SIZE; s2++) if (ht[s2] != EMPTY) occ++;
                    occupied += occ;
                    if (occ == (uint64_t)BIN_SIZE) full_bins++;
                }
            std::cout << "  --- accuracy census ---\n";
            std::cout << "  dropped inserts (bin full): " << g_drops.load() << "\n";
            std::cout << "  full bins / total bins    : " << full_bins << " / " << (N_CU*NBINS) << "\n";
            std::cout << "  occupied slots / capacity : " << occupied << " / " << ((uint64_t)N_CU*NBINS*BIN_SIZE)
                      << "  (load " << 100.0*occupied/((double)N_CU*NBINS*BIN_SIZE) << "%)\n";
        }

        if (print_hist && run == n_runs - 1) {
            std::cout << "\n--- frequency histogram ---\n";
            for (int i = 1; i < 256; i++)
                if (hist[i]) std::cout << i << "\t" << hist[i] << "\n";
        }

        run_e2e_v.push_back(t_run);
    }  /* end run loop */

    if (n_runs > 1) {
        auto med = [](std::vector<double> v) {
            std::sort(v.begin(), v.end());
            int n = v.size();
            return n % 2 == 0 ? (v[n/2-1]+v[n/2])/2.0 : v[n/2];
        };
        std::cout << "\n[stream] Median e2e over " << n_runs << " runs: "
                  << med(run_e2e_v) << " s\n";
    }

    /* ── exp72: phase report ─────────────────────────────────────────── */
    double in_main = elapsed(T_prog, std::chrono::high_resolution_clock::now());
    double sum = 0;
    std::cout << "\n=== SETUP / TEARDOWN PHASE BREAKDOWN (outside batch-loop e2e) ===\n";
    for (auto& pr : g_phase) { sum += pr.second;
        printf("  %-42s %8.3f s\n", pr.first.c_str(), pr.second); }
    printf("  %-42s %8.3f s\n", "[sum of phases above]", sum);
    printf("  %-42s %8.3f s\n", "elapsed inside main() at report time", in_main);
    printf("  %-42s %8llu entries (stride %llu, margin %.2fx)\n",
           "max observed n per (inst,cu)",
           (unsigned long long)g_max_n.load(),
           (unsigned long long)HOST_STRIDE,
           g_max_n.load() ? (double)HOST_STRIDE / (double)g_max_n.load() : 0.0);
    printf("  (process WALL - above) = exit/teardown: free host bufs + device close\n");
    printf("  (host flat_buf+scat_buf total is printed as HOST_STRIDE above; exp73\n"
           "   cut it from 19.3 GB (FE_MAX_PER_CU stride) to ~6 GB at -b 40)\n");
    return 0;
}
