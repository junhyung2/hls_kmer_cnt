#include "kmer_count_opt1.h"

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

/* On-device HT init: write EMPTY to all T1 slots (sequential burst writes) */
static void init_t1(ht_entry_opt1_t* ht_t1)
{
INIT_T1:
    for (uint64_t i = 0; i < T_SIZE_OPT1; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=33554432 avg=2097152
        ht_entry_opt1_t e;
        e.key   = EMPTY_KEY_OPT1;
        e.count = 0;
        e.pad   = 0;
        ht_t1[i] = e;
    }
}

/* On-device HT init: write EMPTY to all T2 slots */
static void init_t2(ht_entry_opt1_t* ht_t2)
{
INIT_T2:
    for (uint64_t i = 0; i < T_SIZE_OPT1; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=33554432 avg=2097152
        ht_entry_opt1_t e;
        e.key   = EMPTY_KEY_OPT1;
        e.count = 0;
        e.pad   = 0;
        ht_t2[i] = e;
    }
}

/* Initialize T1 and T2 simultaneously via DATAFLOW */
static void ht_init(ht_entry_opt1_t* ht_t1, ht_entry_opt1_t* ht_t2)
{
#pragma HLS DATAFLOW
    init_t1(ht_t1);
    init_t2(ht_t2);
}

/* Cuckoo insert: consumes pre-encoded k-mer hash array.
 * Same 2-table direct probe as kmer_count_mlp; PIPELINE II=1 via
 * two independent AXI masters (gmem_t1, gmem_t2).
 * No NO_KMER sentinel check needed — array contains only valid hashes. */
static void cuckoo_insert(
    const uint64_t*   kmers,
    uint64_t          n_kmers,
    ht_entry_opt1_t*  ht_t1,
    ht_entry_opt1_t*  ht_t2,
    uint64_t*         n_unique
) {
    uint64_t unique = 0;

INSERT_LOOP:
    for (uint64_t i = 0; i < n_kmers; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=60000000 avg=30000000
#pragma HLS DEPENDENCE variable=ht_t1 inter false
#pragma HLS DEPENDENCE variable=ht_t2 inter false
        uint64_t h  = kmers[i];
        uint64_t h1 = h & T_MASK_OPT1;
        uint64_t h2 = hash64_c2(h, T_MASK_OPT1);

        ht_entry_opt1_t e1 = ht_t1[h1];
        ht_entry_opt1_t e2 = ht_t2[h2];

        if (e1.key == h) {
            e1.count++; ht_t1[h1] = e1;
        } else if (e2.key == h) {
            e2.count++; ht_t2[h2] = e2;
        } else if (e1.key == EMPTY_KEY_OPT1) {
            ht_entry_opt1_t ne; ne.key = h; ne.count = 1; ne.pad = 0;
            ht_t1[h1] = ne; unique++;
        } else if (e2.key == EMPTY_KEY_OPT1) {
            ht_entry_opt1_t ne; ne.key = h; ne.count = 1; ne.pad = 0;
            ht_t2[h2] = ne; unique++;
        }
        /* else: double collision → drop (P < 0.09% at HW load factor) */
    }
    *n_unique = unique;
}

extern "C" void kmer_count_opt1(
    const uint64_t*   kmers,
    ht_entry_opt1_t*  ht_t1,
    ht_entry_opt1_t*  ht_t2,
    uint64_t          n_kmers,
    uint64_t*         n_unique
) {
#pragma HLS INTERFACE m_axi port=kmers   offset=slave bundle=gmem0 \
        max_read_burst_length=256 num_read_outstanding=16
#pragma HLS INTERFACE m_axi port=ht_t1  offset=slave bundle=gmem_t1 \
        max_read_burst_length=2 max_write_burst_length=256 \
        num_read_outstanding=32 num_write_outstanding=32
#pragma HLS INTERFACE m_axi port=ht_t2  offset=slave bundle=gmem_t2 \
        max_read_burst_length=2 max_write_burst_length=256 \
        num_read_outstanding=32 num_write_outstanding=32
#pragma HLS INTERFACE m_axi port=n_unique offset=slave bundle=gmem3

#pragma HLS INTERFACE s_axilite port=kmers    bundle=control
#pragma HLS INTERFACE s_axilite port=ht_t1    bundle=control
#pragma HLS INTERFACE s_axilite port=ht_t2    bundle=control
#pragma HLS INTERFACE s_axilite port=n_kmers  bundle=control
#pragma HLS INTERFACE s_axilite port=n_unique bundle=control
#pragma HLS INTERFACE s_axilite port=return   bundle=control

    /* Step 1: initialize HT on-device (T1 and T2 in parallel via DATAFLOW) */
    ht_init(ht_t1, ht_t2);

    /* Step 2: insert pre-encoded k-mers */
    cuckoo_insert(kmers, n_kmers, ht_t1, ht_t2, n_unique);
}
