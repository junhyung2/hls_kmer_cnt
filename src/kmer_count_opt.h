#ifndef KMER_COUNT_OPT_H
#define KMER_COUNT_OPT_H

#include <stdint.h>
#include <hls_stream.h>
#include <ap_int.h>

/* Per-CU hash table sizing.
 * With 4 CUs routing by bottom 2 bits of hash, each CU handles ~1/4 of all
 * unique k-mers.  Table is split: T1 = [0, T_SIZE), T2 = [T_SIZE, HT_SIZE).
 *
 *   SW_EMU / HW_EMU : 2^22 = 4M entries x 16B =  64MB per CU (fast sim)
 *   HW              : 2^26 = 64M entries x 16B = 1GB  per CU (4 CUs = 4GB)
 */
#if defined(SW_EMU) || defined(HW_EMU)
  #define HT_BITS_OPT  22
#else
  #define HT_BITS_OPT  26
#endif

#define HT_SIZE_OPT  (1ULL << HT_BITS_OPT)
#define HT_MASK_OPT  (HT_SIZE_OPT - 1ULL)

/* Cuckoo split: T1 = lower half, T2 = upper half */
#define T_SIZE_OPT   (HT_SIZE_OPT / 2ULL)
#define T_MASK_OPT   (T_SIZE_OPT  - 1ULL)

#define EMPTY_KEY_OPT  0xFFFFFFFFFFFFFFFFULL
#define NO_KMER_OPT    0xFFFFFFFFFFFFFFFFULL

/* Max linear-probe steps in T1 fallback (rare: both T1 and T2 occupied) */
#define MAX_PROBES_OPT  128

/* Number of CUs */
#define N_CU  4

/* Stream depth: SW_EMU DATAFLOW is sequential → buffer full parse output */
#if defined(SW_EMU)
  #define STREAM_DEPTH_OPT  5200000
#else
  #define STREAM_DEPTH_OPT  4096
#endif

typedef struct {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
} ht_entry_opt_t;

/* Top-level kernel (instantiated N_CU times):
 *   seq1, seq2  read-only input (shared across all CU instances)
 *   ht          per-CU hash table (separate buffer per CU)
 *   n_unique    per-CU output
 *   cu_id       [0,3] selects which quarter of k-mers this CU processes
 */
extern "C" void kmer_count_opt(
    const uint8_t*   seq1,
    const uint8_t*   seq2,
    ht_entry_opt_t*  ht,
    uint64_t         len1,
    uint64_t         len2,
    int              k,
    uint64_t*        n_unique,
    int              cu_id
);

#endif /* KMER_COUNT_OPT_H */
