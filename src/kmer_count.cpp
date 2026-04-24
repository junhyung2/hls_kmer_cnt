#include "kmer_count.h"
#include <string.h>

/* NT4 lookup: A=0 C=1 G=2 T/U=3 other=4
 * 256-entry ROM -> inferred as LUTROM by HLS */
static const uint8_t NT4[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  /* A=0 C=1 G=2 */
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,  /* T=3 U=3 */
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  /* a=0 c=1 g=2 */
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,  /* t=3 u=3 */
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4
};

/* kc-c4 invertible hash: uniform distribution over [0, mask] */
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

/* Stage 1: FASTQ parser + canonical k-mer extractor
 *
 * Reads both FASTQ files sequentially in a single function.
 * (DATAFLOW rule: each m_axi input must be accessed by exactly one process.)
 *
 * Emits exactly (len1 + len2) values into 'out':
 *   valid k-mer  ->  hash64(canonical, km_mask)
 *   other byte   ->  NO_KMER  (0xFFFF...FF)
 *
 * 1-to-1 byte-to-stream mapping gives hash_insert a fixed loop bound.
 *
 * FASTQ record_line:  0=@header  1=sequence  2=+line  3=quality */
static void parse_extract_both(
    const uint8_t* seq1,
    uint64_t       len1,
    const uint8_t* seq2,
    uint64_t       len2,
    int            k,
    uint64_t       km_mask,
    int            shift,
    hls::stream<uint64_t>& out
)
{
    uint64_t x0    = 0;
    uint64_t x1    = 0;
    int      l     = 0;
    int      rline = 0;

SEQ1_LOOP:
    for (uint64_t i = 0; i < len1; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=209715200 avg=2621440
        uint8_t  raw  = seq1[i];
        uint64_t emit = NO_KMER;

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
                    emit = hash64(y, km_mask);
                }
            } else {
                l = 0; x0 = 0; x1 = 0;
            }
        }
        out.write(emit);
    }

    /* Reset sliding window state between files */
    x0 = 0; x1 = 0; l = 0; rline = 0;

SEQ2_LOOP:
    for (uint64_t i = 0; i < len2; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=209715200 avg=2621440
        uint8_t  raw  = seq2[i];
        uint64_t emit = NO_KMER;

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
                    emit = hash64(y, km_mask);
                }
            } else {
                l = 0; x0 = 0; x1 = 0;
            }
        }
        out.write(emit);
    }
}

/* Stage 2: stream consumer -> hash table updater
 *
 * Reads exactly 'total_len' values from 'in'.
 * Skips NO_KMER entries; inserts/increments valid k-mers in the hash table.
 * Open addressing with linear probing, bounded at MAX_PROBES steps. */
static void hash_insert(
    hls::stream<uint64_t>& in,
    uint64_t    total_len,
    ht_entry_t* ht,
    uint64_t*   n_unique
)
{
    uint64_t unique = 0;

OUTER_LOOP:
    for (uint64_t i = 0; i < total_len; i++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=419430400 avg=5242880
        uint64_t h = in.read();
        if (h == NO_KMER) continue;

        uint64_t bucket = h & HT_MASK;
        bool     done   = false;

    PROBE_LOOP:
        for (int p = 0; p < MAX_PROBES; p++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=2
            if (done) break;
            ht_entry_t e = ht[bucket];

            if (e.key == EMPTY_KEY) {
                ht_entry_t ne;
                ne.key   = h;
                ne.count = 1;
                ne.pad   = 0;
                ht[bucket] = ne;
                unique++;
                done = true;
            } else if (e.key == h) {
                e.count++;
                ht[bucket] = e;
                done = true;
            } else {
                bucket = (bucket + 1) & HT_MASK;
            }
        }
    }

    *n_unique = unique;
}

/* Top-level HLS kernel */
extern "C" void kmer_count(
    const uint8_t* seq1,
    const uint8_t* seq2,
    ht_entry_t*    ht,
    uint64_t       len1,
    uint64_t       len2,
    int            k,
    uint64_t*      n_unique
)
{
/* seq1/seq2: burst reads for sequential FASTQ scan */
#pragma HLS INTERFACE m_axi port=seq1     offset=slave bundle=gmem0 \
        max_read_burst_length=256 num_read_outstanding=16
#pragma HLS INTERFACE m_axi port=seq2     offset=slave bundle=gmem1 \
        max_read_burst_length=256 num_read_outstanding=16
/* ht: random RMW; burst=2 reflects random access pattern */
#pragma HLS INTERFACE m_axi port=ht       offset=slave bundle=gmem2 \
        max_read_burst_length=2   max_write_burst_length=2 \
        num_read_outstanding=32   num_write_outstanding=32
/* n_unique: single 64-bit write at end */
#pragma HLS INTERFACE m_axi port=n_unique offset=slave bundle=gmem3

#pragma HLS INTERFACE s_axilite port=seq1     bundle=control
#pragma HLS INTERFACE s_axilite port=seq2     bundle=control
#pragma HLS INTERFACE s_axilite port=ht       bundle=control
#pragma HLS INTERFACE s_axilite port=len1     bundle=control
#pragma HLS INTERFACE s_axilite port=len2     bundle=control
#pragma HLS INTERFACE s_axilite port=k        bundle=control
#pragma HLS INTERFACE s_axilite port=n_unique bundle=control
#pragma HLS INTERFACE s_axilite port=return   bundle=control

    uint64_t km_mask = (1ULL << (2 * k)) - 1;
    int      shift   = (k - 1) * 2;
    uint64_t total   = len1 + len2;

    /* Inter-stage stream.
     * Depth sized per target (see STREAM_DEPTH in kmer_count.h):
     *   SW_EMU  : 5,200,000 (covers full output in sequential simulation)
     *   HW/EMU  : 4,096     (back-pressure buffer for concurrent execution) */
    hls::stream<uint64_t> kmer_stream;
#if defined(SW_EMU)
#pragma HLS STREAM variable=kmer_stream depth=5200000
#else
#pragma HLS STREAM variable=kmer_stream depth=4096
#endif

/* DATAFLOW: parse_extract_both and hash_insert run concurrently in HW.
 * In SW_EMU they execute sequentially; stream depth covers full output. */
#pragma HLS DATAFLOW

    parse_extract_both(seq1, len1, seq2, len2, k, km_mask, shift, kmer_stream);
    hash_insert(kmer_stream, total, ht, n_unique);
}
