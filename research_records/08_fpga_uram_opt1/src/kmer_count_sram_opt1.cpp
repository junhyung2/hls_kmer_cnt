#include "kmer_count_sram_opt1.h"

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

extern "C" void kmer_count_sram_opt1(
    const uint64_t*  kmers,
    const uint64_t*  seg_ends,
    uint64_t         n_kmers,
    uint64_t*        n_unique
) {
#pragma HLS INTERFACE m_axi port=kmers    offset=slave bundle=gmem0 \
        max_read_burst_length=256 num_read_outstanding=16
#pragma HLS INTERFACE m_axi port=seg_ends offset=slave bundle=gmem1 \
        max_read_burst_length=16  num_read_outstanding=1
#pragma HLS INTERFACE m_axi port=n_unique offset=slave bundle=gmem2

#pragma HLS INTERFACE s_axilite port=kmers    bundle=control
#pragma HLS INTERFACE s_axilite port=seg_ends bundle=control
#pragma HLS INTERFACE s_axilite port=n_kmers  bundle=control
#pragma HLS INTERFACE s_axilite port=n_unique bundle=control
#pragma HLS INTERFACE s_axilite port=return   bundle=control

    ht_entry_sram_opt1_t ht_t1[SRAM_T_SIZE];
    ht_entry_sram_opt1_t ht_t2[SRAM_T_SIZE];
#pragma HLS BIND_STORAGE variable=ht_t1 type=RAM_2P impl=URAM
#pragma HLS BIND_STORAGE variable=ht_t2 type=RAM_2P impl=URAM

    /* ---------------------------------------------------------------
     * One-time initialization: mark every slot as belonging to no
     * segment (EMPTY_SEG_ID = 0xFFFF).  Runs once in 2×131072 cycles
     * (~0.87 ms at 300 MHz) instead of the per-segment clear that
     * cost 2×131072×256 ≈ 67M cycles in the 07 design.
     * --------------------------------------------------------------- */
    ht_entry_sram_opt1_t empty_e;
    empty_e.key    = EMPTY_KEY_SRAM;
    empty_e.count  = 0;
    empty_e.seg_id = EMPTY_SEG_ID;
    empty_e.pad    = 0;

INIT_T1:
    for (int i = 0; i < SRAM_T_SIZE; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=131072 max=131072 avg=131072
        ht_t1[i] = empty_e;
    }
INIT_T2:
    for (int i = 0; i < SRAM_T_SIZE; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=131072 max=131072 avg=131072
        ht_t2[i] = empty_e;
    }

    /* Cache seg_ends in registers */
    uint64_t se[N_SEG];
#pragma HLS ARRAY_PARTITION variable=se complete
READ_SE:
    for (int s = 0; s < N_SEG; s++) {
#pragma HLS PIPELINE II=1
        se[s] = seg_ends[s];
    }

    uint64_t total_unique = 0;

    /* ---------------------------------------------------------------
     * Segment loop — NO clear between segments.
     *
     * Lazy-clear: each slot carries seg_id (the segment that last
     * wrote it).  A slot is "valid for current segment s" iff
     * e.seg_id == s.  Otherwise it is stale and treated as empty.
     *
     * Correctness: segments are processed in strict order 0..255.
     * A slot written by segment s' < s has seg_id = s' ≠ s → stale.
     * The one-time INIT guarantees seg_id = 0xFFFF ≠ any s ∈ [0,255]
     * before the first segment, so segment 0 sees a completely empty
     * table.
     *
     * Timing (HW, 300 MHz):
     *   init     : 2×131072 / 300M ≈ 0.87 ms  (one-time)
     *   insert   : ~39K k-mer/seg/CU × II=1 / 300M ≈ 0.13 ms/seg
     *   × 256 seg ≈ 33 ms total  (vs 07: 290 ms — ~8.8× speedup)
     * --------------------------------------------------------------- */
SEG_LOOP:
    for (int seg = 0; seg < N_SEG; seg++) {
#pragma HLS LOOP_TRIPCOUNT min=256 max=256 avg=256

        uint64_t seg_start = (seg == 0) ? 0ULL : se[seg - 1];
        uint64_t seg_end   = se[seg];
        uint64_t seg_unique = 0;
        uint16_t cur_seg    = (uint16_t)seg;

INSERT_LOOP:
        for (uint64_t j = seg_start; j < seg_end; j++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=2000000 avg=39000
#pragma HLS DEPENDENCE variable=ht_t1 inter false
#pragma HLS DEPENDENCE variable=ht_t2 inter false

            uint64_t h  = kmers[j];
            uint64_t h1 = (h >> N_ROUTE_BITS) & SRAM_T_MASK;
            uint64_t h2 = hash64_c2(h, SRAM_T_MASK);

            ht_entry_sram_opt1_t e1 = ht_t1[h1];
            ht_entry_sram_opt1_t e2 = ht_t2[h2];

            bool t1_mine = (e1.seg_id == cur_seg);
            bool t2_mine = (e2.seg_id == cur_seg);

            if (t1_mine && e1.key == h) {
                e1.count++;
                ht_t1[h1] = e1;
            } else if (t2_mine && e2.key == h) {
                e2.count++;
                ht_t2[h2] = e2;
            } else if (!t1_mine) {
                /* T1 slot is stale — overwrite */
                ht_entry_sram_opt1_t ne;
                ne.key = h; ne.count = 1; ne.seg_id = cur_seg; ne.pad = 0;
                ht_t1[h1] = ne;
                seg_unique++;
            } else if (!t2_mine) {
                /* T2 slot is stale — overwrite */
                ht_entry_sram_opt1_t ne;
                ne.key = h; ne.count = 1; ne.seg_id = cur_seg; ne.pad = 0;
                ht_t2[h2] = ne;
                seg_unique++;
            }
            /* else: both T1 and T2 occupied by cur_seg with diff key → drop */
        }

        total_unique += seg_unique;
    }

    *n_unique = total_unique;
}
