#pragma once
#include <stdint.h>
#include <hls_stream.h>
#include <ap_int.h>

/* ---------------------------------------------------------------
 * Routing bits:
 *   bits [N_CU_BITS-1 : 0]                        → CU  selection (4 CUs)
 *   bits [N_CU_BITS+N_SEG_BITS-1 : N_CU_BITS]     → seg selection (256 segs)
 *   → 10 routing bits consumed before URAM hash index
 * --------------------------------------------------------------- */
#define N_CU_BITS      2                         /* 4 CUs per card          */
#define N_SEG_BITS     8                         /* 256 segments per CU     */
#define N_ROUTE_BITS   (N_CU_BITS + N_SEG_BITS) /* 10 bits total           */
#define N_SEG          (1 << N_SEG_BITS)         /* 256                     */
#define N_CU_SRAM      4

/* ---------------------------------------------------------------
 * URAM hash table sizing
 *   HW  : HT_BITS=18 → T_SIZE=131072 entries per T per segment
 *          Vivado measured: 72 URAMs per CU (T1+T2)
 *          Total 4 CU: 288 / 640 URAMs (45%)
 *   EMU : HT_BITS=14 → T_SIZE=8192  (fits in BRAM for fast emulation)
 *
 *   Accuracy analysis (HW, 1 card 4 CU, N_SEG=256):
 *     unique/seg/CU = 9M / (4 × 256)  ≈ 8,800
 *     T1+T2 slots   = 2 × 131072      = 262,144
 *     drop_fraction ≈ N²/(12×T²)      ≈ 0.038%
 *     → comparable accuracy to opt1 (0.05% drop)
 *
 *   Previous N_SEG=16 caused 7% drop (load factor 0.54 → too high).
 *   N_SEG=256 reduces per-segment unique to 8.8K → drop negligible.
 * --------------------------------------------------------------- */
#if defined(SW_EMU) || defined(HW_EMU) || defined(CSIM)
  #define SRAM_HT_BITS  14
#else
  #define SRAM_HT_BITS  18
#endif

#define SRAM_T_SIZE  (1 << (SRAM_HT_BITS - 1))   /* entries per T (131072 HW) */
#define SRAM_T_MASK  (SRAM_T_SIZE - 1ULL)

#define EMPTY_KEY_SRAM  0xFFFFFFFFFFFFFFFFULL

struct ht_entry_sram_t {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
};

extern "C" void kmer_count_sram(
    const uint64_t*       kmers,    /* segments 0..N_SEG-1 concatenated      */
    const uint64_t*       seg_ends, /* [N_SEG] cumulative end offsets         */
    uint64_t              n_kmers,  /* total k-mers for this CU               */
    uint64_t*             n_unique  /* output: total unique k-mers            */
);
