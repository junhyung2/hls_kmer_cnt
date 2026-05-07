#pragma once
#include <stdint.h>
#include <hls_stream.h>
#include <ap_int.h>

/* ---------------------------------------------------------------
 * Routing layout (bits of canonical-kmer hash h):
 *   [1:0]   cu_id   — 4 CUs per card
 *   [9:2]   seg_id  — 256 segments per CU
 *   [26:10] T1 idx  — 17-bit T1 index (131072 slots)
 *   [41:27] card    — (h>>27) % 3  → avoids overlap with T1 bits
 * --------------------------------------------------------------- */
#define N_CU_BITS      2
#define N_SEG_BITS     8
#define N_ROUTE_BITS   (N_CU_BITS + N_SEG_BITS)   /* 10 */
#define N_SEG          (1 << N_SEG_BITS)           /* 256 */
#define N_CU_SRAM      4

#if defined(SW_EMU) || defined(HW_EMU) || defined(CSIM)
  #define SRAM_HT_BITS  14
#else
  #define SRAM_HT_BITS  18
#endif

#define SRAM_T_SIZE  (1 << (SRAM_HT_BITS - 1))    /* 131072 HW */
#define SRAM_T_MASK  (SRAM_T_SIZE - 1ULL)

/* card routing shift: T1 uses bits[26:10], card uses bits[41:27] */
#define CARD_ROUTE_SHIFT  (N_ROUTE_BITS + SRAM_HT_BITS - 1)  /* 27 HW */

#define EMPTY_KEY_SRAM  0xFFFFFFFFFFFFFFFFULL
#define EMPTY_SEG_ID    0xFFFFU               /* marks uninitialized slot */

/* seg_id replaces pad[31:16]: struct size unchanged at 16 B */
struct ht_entry_sram_opt1_t {
    uint64_t key;
    uint32_t count;
    uint16_t seg_id;
    uint16_t pad;
};

extern "C" void kmer_count_sram_opt1(
    const uint64_t*  kmers,
    const uint64_t*  seg_ends,
    uint64_t         n_kmers,
    uint64_t*        n_unique
);
