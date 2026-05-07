#ifndef KMER_COUNT_MLP_H
#define KMER_COUNT_MLP_H

#include <stdint.h>
#include <hls_stream.h>
#include <ap_int.h>

/* MLP-optimized kmer_count: T1 and T2 are SEPARATE AXI master bundles
 * (bundle=gmem_t1, bundle=gmem_t2).  This breaks the T1→T2 sequential
 * dependency inside cuckoo_insert, enabling #pragma HLS PIPELINE II=1
 * on the outer loop.  With num_read_outstanding=32, both masters sustain
 * 32 in-flight requests each → 64 total in-flight → full HBM latency hidden.
 *
 * Key difference from kmer_count_opt (single ht bundle):
 *   Old: e1=ht[h1] → wait → e2=ht[h2] → wait → write  (serial, II~60)
 *   New: e1=ht_t1[h1] || e2=ht_t2[h2] issued same cycle  (parallel, II=1)
 *
 * Hash table sizing (each T array holds T_SIZE_MLP entries):
 *   CSIM / SW_EMU / HW_EMU : HT_BITS=22 → T_SIZE=2^21=2M entries = 32 MB per T
 *   HW                     : HT_BITS=26 → T_SIZE=2^25=32M entries = 512 MB per T
 *
 * HBM layout (4 CUs, HW target):
 *   seq1     → HBM[0]         (shared read-only)
 *   seq2     → HBM[1]         (shared read-only)
 *   CU1.ht_t1 → HBM[2:3]      CU1.ht_t2 → HBM[4:5]
 *   CU2.ht_t1 → HBM[6:7]      CU2.ht_t2 → HBM[8:9]
 *   CU3.ht_t1 → HBM[10:11]    CU3.ht_t2 → HBM[12:13]
 *   CU4.ht_t1 → HBM[14:15]    CU4.ht_t2 → HBM[16:17]
 *   n_unique  → HBM[18]        (shared write)
 *
 * HMSS port count: 4+4+4+4+4 = 20 (well within U50 limit of 32)
 */

#if defined(SW_EMU) || defined(HW_EMU) || defined(CSIM)
  #define HT_BITS_MLP  22
#else
  #define HT_BITS_MLP  26
#endif

/* T1 and T2 are separate arrays, each of size T_SIZE_MLP */
#define T_SIZE_MLP   (1ULL << (HT_BITS_MLP - 1))
#define T_MASK_MLP   (T_SIZE_MLP - 1ULL)

#define EMPTY_KEY_MLP  0xFFFFFFFFFFFFFFFFULL
#define NO_KMER_MLP    0xFFFFFFFFFFFFFFFFULL

#define N_CU_MLP  4

typedef struct {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
} ht_entry_mlp_t;

/* Top-level kernel (instantiated N_CU_MLP times):
 *   seq1, seq2  shared read-only input
 *   ht_t1       per-CU T1 region (AXI master gmem_t1)
 *   ht_t2       per-CU T2 region (AXI master gmem_t2, separate from ht_t1)
 *   n_unique    per-CU unique k-mer count output
 *   cu_id       [0,3]: selects which quarter of k-mers this CU processes
 */
extern "C" void kmer_count_mlp(
    const uint8_t*   seq1,
    const uint8_t*   seq2,
    ht_entry_mlp_t*  ht_t1,
    ht_entry_mlp_t*  ht_t2,
    uint64_t         len1,
    uint64_t         len2,
    int              k,
    uint64_t*        n_unique,
    int              cu_id
);

#endif /* KMER_COUNT_MLP_H */
