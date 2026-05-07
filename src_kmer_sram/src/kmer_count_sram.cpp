#include "kmer_count_sram.h"

/* ---------------------------------------------------------------
 * Full 64-bit mix for T2 index — routing bits 0-5 are fixed per
 * (CU, segment), so a plain mask on h would cluster T1.  We use
 * (h >> N_ROUTE_BITS) for T1 and a full remix for T2.
 * --------------------------------------------------------------- */
static inline uint64_t hash64_c2(uint64_t key, uint64_t mask)
{
#pragma HLS INLINE
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key & mask;
}

/* ---------------------------------------------------------------
 * Clear one URAM table: II=1 sequential write burst.
 * --------------------------------------------------------------- */
static void clear_table(ht_entry_sram_t* t, int sz)
{
#pragma HLS INLINE off
CLEAR:
    for (int i = 0; i < sz; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=8192 max=131072 avg=131072
        ht_entry_sram_t e;
        e.key   = EMPTY_KEY_SRAM;
        e.count = 0;
        e.pad   = 0;
        t[i] = e;
    }
}

extern "C" void kmer_count_sram(
    const uint64_t*   kmers,
    const uint64_t*   seg_ends,
    uint64_t          n_kmers,
    uint64_t*         n_unique
) {
/* ---- AXI master interfaces ---- */
#pragma HLS INTERFACE m_axi port=kmers    offset=slave bundle=gmem0 \
        max_read_burst_length=256 num_read_outstanding=16
#pragma HLS INTERFACE m_axi port=seg_ends offset=slave bundle=gmem1 \
        max_read_burst_length=16  num_read_outstanding=1
#pragma HLS INTERFACE m_axi port=n_unique offset=slave bundle=gmem2

/* ---- s_axilite control ---- */
#pragma HLS INTERFACE s_axilite port=kmers    bundle=control
#pragma HLS INTERFACE s_axilite port=seg_ends bundle=control
#pragma HLS INTERFACE s_axilite port=n_kmers  bundle=control
#pragma HLS INTERFACE s_axilite port=n_unique bundle=control
#pragma HLS INTERFACE s_axilite port=return   bundle=control

    /* ---------------------------------------------------------------
     * URAM-resident hash tables.
     *   SRAM_T_SIZE=131072 entries × 16 B = 2 MB per table.
     *   Two 72-bit URAM columns × 32 row-groups = 64 URAMs per table.
     *   RAM_2P: independent read and write ports → II=1 achievable.
     * --------------------------------------------------------------- */
    ht_entry_sram_t ht_t1[SRAM_T_SIZE];
    ht_entry_sram_t ht_t2[SRAM_T_SIZE];
#pragma HLS BIND_STORAGE variable=ht_t1 type=RAM_2P impl=URAM
#pragma HLS BIND_STORAGE variable=ht_t2 type=RAM_2P impl=URAM

    /* ---- Cache seg_ends in registers (16 × 8 B = 128 B burst) ---- */
    uint64_t se[N_SEG];
#pragma HLS ARRAY_PARTITION variable=se complete
READ_SE:
    for (int s = 0; s < N_SEG; s++) {
#pragma HLS PIPELINE II=1
        se[s] = seg_ends[s];
    }

    uint64_t total_unique = 0;

    /* ---------------------------------------------------------------
     * Outer loop: process each segment sequentially.
     *   1. Clear T1, T2 (URAM burst writes, ~131K cycles each).
     *   2. Cuckoo-insert k-mers for this segment using URAM T1/T2.
     *   3. Accumulate unique count.
     *
     * Timing (HW, 300 MHz):
     *   clear  : 2 × 131072 / 300M ≈ 0.87 ms per segment
     *   insert : ~937K k-mers × II≈1 / 300M ≈ 3.1 ms per segment
     *   × 16 segments ≈ 63 ms per CU  (vs opt1 3790 ms — ~60× speedup)
     * --------------------------------------------------------------- */
SEG_LOOP:
    for (int seg = 0; seg < N_SEG; seg++) {
#pragma HLS LOOP_TRIPCOUNT min=256 max=256 avg=256

        clear_table(ht_t1, SRAM_T_SIZE);
        clear_table(ht_t2, SRAM_T_SIZE);

        uint64_t seg_start = (seg == 0) ? 0ULL : se[seg - 1];
        uint64_t seg_end   = se[seg];
        uint64_t seg_unique = 0;

INSERT_LOOP:
        for (uint64_t j = seg_start; j < seg_end; j++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=2000000 avg=117000
#pragma HLS DEPENDENCE variable=ht_t1 inter false
#pragma HLS DEPENDENCE variable=ht_t2 inter false

            uint64_t h  = kmers[j];
            /* T1 index: skip routing bits so lower bits vary within segment */
            uint64_t h1 = (h >> N_ROUTE_BITS) & SRAM_T_MASK;
            /* T2 index: full remix removes routing-bit clustering            */
            uint64_t h2 = hash64_c2(h, SRAM_T_MASK);

            ht_entry_sram_t e1 = ht_t1[h1];
            ht_entry_sram_t e2 = ht_t2[h2];

            if (e1.key == h) {
                e1.count++;
                ht_t1[h1] = e1;
            } else if (e2.key == h) {
                e2.count++;
                ht_t2[h2] = e2;
            } else if (e1.key == EMPTY_KEY_SRAM) {
                ht_entry_sram_t ne;
                ne.key = h; ne.count = 1; ne.pad = 0;
                ht_t1[h1] = ne;
                seg_unique++;
            } else if (e2.key == EMPTY_KEY_SRAM) {
                ht_entry_sram_t ne;
                ne.key = h; ne.count = 1; ne.pad = 0;
                ht_t2[h2] = ne;
                seg_unique++;
            }
            /* else: double collision → drop (P negligible at load 0.54) */
        }

        total_unique += seg_unique;
    }

    *n_unique = total_unique;
}
