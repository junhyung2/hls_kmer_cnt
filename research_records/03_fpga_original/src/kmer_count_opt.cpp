#include "kmer_count_opt.h"
#include <string.h>

/* NT4: A=0 C=1 G=2 T/U=3 other=4 (256-entry ROM) */
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

/* Primary hash: kc-c4 invertible hash (same as original kmer_count) */
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
#pragma HLS INLINE
    key = (~key + (key << 21)) & mask;
    key =   key ^ (key >> 24);
    key = ((key + (key <<  3)) + (key <<  8)) & mask;
    key =   key ^ (key >> 14);
    key = ((key + (key <<  2)) + (key <<  4)) & mask;
    key =   key ^ (key >> 28);
    key =  (key + (key << 31)) & mask;
    return key;
}

/* Secondary hash for cuckoo T2 slot (Murmur3 64-bit finalizer) */
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

/* Stage 1: FASTQ parser + canonical k-mer extractor with CU routing.
 *
 * Scans both FASTQ files sequentially.  For each valid k-mer hash h:
 *   if (h & 3) == cu_id  →  emit h (this CU's k-mer)
 *   else                 →  emit NO_KMER_OPT (skip)
 *
 * Always emits exactly (len1 + len2) values so cuckoo_insert has a fixed
 * trip count and DATAFLOW back-pressure is bounded. */
static void parse_extract_cu(
    const uint8_t*       seq1,
    uint64_t             len1,
    const uint8_t*       seq2,
    uint64_t             len2,
    int                  k,
    uint64_t             km_mask,
    int                  shift,
    int                  cu_id,
    hls::stream<uint64_t>& out
) {
    uint64_t x0    = 0, x1 = 0;
    int      l     = 0, rline = 0;

SEQ1_LOOP:
    for (uint64_t i = 0; i < len1; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=209715200 avg=2621440
        uint8_t  raw  = seq1[i];
        uint64_t emit = NO_KMER_OPT;

        if (raw == (uint8_t)'\n') {
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
                    if ((int)(h & 3) == cu_id) emit = h;
                }
            } else { l = 0; x0 = 0; x1 = 0; }
        }
        out.write(emit);
    }

    x0 = 0; x1 = 0; l = 0; rline = 0;

SEQ2_LOOP:
    for (uint64_t i = 0; i < len2; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=209715200 avg=2621440
        uint8_t  raw  = seq2[i];
        uint64_t emit = NO_KMER_OPT;

        if (raw == (uint8_t)'\n') {
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
                    if ((int)(h & 3) == cu_id) emit = h;
                }
            } else { l = 0; x0 = 0; x1 = 0; }
        }
        out.write(emit);
    }
}

/* Stage 2: Cuckoo hash insert.
 *
 * Two-choice hashing with linear-probe fallback:
 *   T1 slot h1 = h & T_MASK_OPT           → ht[h1]
 *   T2 slot h2 = T_SIZE_OPT + c2(h)       → ht[h2]
 *
 * Lookup/insert order:
 *   1. Check T1[h1]: match → increment; else
 *   2. Check T2[h2]: match → increment; else
 *   3. T1[h1] empty → insert in T1; else
 *   4. T2[h2] empty → insert in T2; else
 *   5. Linear probe T1 only (wrap within [0, T_SIZE_OPT))
 *
 * Correctness: keys inserted by linear probe are in T1 region only.
 * Lookup correctly finds them because the probe re-starts from h1+1.
 * T2 slots are accessed only at the fixed h2 address, never via linear scan. */
static void cuckoo_insert(
    hls::stream<uint64_t>& in,
    uint64_t               total_len,
    ht_entry_opt_t*        ht,
    uint64_t*              n_unique
) {
    uint64_t unique = 0;

OUTER_LOOP:
    for (uint64_t i = 0; i < total_len; i++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=419430400 avg=5242880
        uint64_t h = in.read();
        if (h == NO_KMER_OPT) continue;

        uint64_t h1 = h & T_MASK_OPT;
        uint64_t h2 = T_SIZE_OPT + hash64_c2(h, T_MASK_OPT);

        /* --- Phase 1: T1 direct probe --- */
        ht_entry_opt_t e1 = ht[h1];
        if (e1.key == h) {
            e1.count++;
            ht[h1] = e1;
            continue;
        }

        /* --- Phase 2: T2 direct probe --- */
        ht_entry_opt_t e2 = ht[h2];
        if (e2.key == h) {
            e2.count++;
            ht[h2] = e2;
            continue;
        }

        /* --- Phase 3: insert into first available slot --- */
        if (e1.key == EMPTY_KEY_OPT) {
            ht_entry_opt_t ne; ne.key = h; ne.count = 1; ne.pad = 0;
            ht[h1] = ne;
            unique++;
            continue;
        }
        if (e2.key == EMPTY_KEY_OPT) {
            ht_entry_opt_t ne; ne.key = h; ne.count = 1; ne.pad = 0;
            ht[h2] = ne;
            unique++;
            continue;
        }

        /* --- Phase 4: linear probe in T1 (both slots occupied) --- */
        uint64_t bucket = (h1 + 1) & T_MASK_OPT;
    PROBE_LOOP:
        for (int p = 0; p < MAX_PROBES_OPT; p++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=2
            ht_entry_opt_t ep = ht[bucket];
            if (ep.key == h) {
                ep.count++;
                ht[bucket] = ep;
                break;
            }
            if (ep.key == EMPTY_KEY_OPT) {
                ht_entry_opt_t ne; ne.key = h; ne.count = 1; ne.pad = 0;
                ht[bucket] = ne;
                unique++;
                break;
            }
            bucket = (bucket + 1) & T_MASK_OPT;
        }
    }

    *n_unique = unique;
}

/* Top-level HLS kernel */
extern "C" void kmer_count_opt(
    const uint8_t*   seq1,
    const uint8_t*   seq2,
    ht_entry_opt_t*  ht,
    uint64_t         len1,
    uint64_t         len2,
    int              k,
    uint64_t*        n_unique,
    int              cu_id
) {
/* seq1/seq2: sequential burst reads (all 4 CU instances share same HBM banks) */
#pragma HLS INTERFACE m_axi port=seq1     offset=slave bundle=gmem0 \
        max_read_burst_length=256 num_read_outstanding=16
#pragma HLS INTERFACE m_axi port=seq2     offset=slave bundle=gmem1 \
        max_read_burst_length=256 num_read_outstanding=16
/* ht: random RMW; burst=2 for random access */
#pragma HLS INTERFACE m_axi port=ht       offset=slave bundle=gmem2 \
        max_read_burst_length=2   max_write_burst_length=2 \
        num_read_outstanding=32   num_write_outstanding=32
/* n_unique: single write */
#pragma HLS INTERFACE m_axi port=n_unique offset=slave bundle=gmem3

#pragma HLS INTERFACE s_axilite port=seq1     bundle=control
#pragma HLS INTERFACE s_axilite port=seq2     bundle=control
#pragma HLS INTERFACE s_axilite port=ht       bundle=control
#pragma HLS INTERFACE s_axilite port=len1     bundle=control
#pragma HLS INTERFACE s_axilite port=len2     bundle=control
#pragma HLS INTERFACE s_axilite port=k        bundle=control
#pragma HLS INTERFACE s_axilite port=n_unique bundle=control
#pragma HLS INTERFACE s_axilite port=cu_id    bundle=control
#pragma HLS INTERFACE s_axilite port=return   bundle=control

    uint64_t km_mask = (1ULL << (2 * k)) - 1;
    int      shift   = (k - 1) * 2;
    uint64_t total   = len1 + len2;

    hls::stream<uint64_t> kmer_stream;
#if defined(SW_EMU)
#pragma HLS STREAM variable=kmer_stream depth=5200000
#else
#pragma HLS STREAM variable=kmer_stream depth=4096
#endif

#pragma HLS DATAFLOW
    parse_extract_cu(seq1, len1, seq2, len2, k, km_mask, shift, cu_id, kmer_stream);
    cuckoo_insert(kmer_stream, total, ht, n_unique);
}
