#pragma once
#include <stdint.h>

#ifndef __SYNTHESIS__
#include <hls_stream.h>
#endif
#include <hls_stream.h>
#include <ap_int.h>

/* HT_BITS_OPT1: set per target in Makefile / CSIM define */
#if defined(SW_EMU) || defined(HW_EMU) || defined(CSIM)
  #define HT_BITS_OPT1  22
#else
  #define HT_BITS_OPT1  26
#endif

#define T_SIZE_OPT1   (1ULL << (HT_BITS_OPT1 - 1))
#define T_MASK_OPT1   (T_SIZE_OPT1 - 1ULL)
#define EMPTY_KEY_OPT1  0xFFFFFFFFFFFFFFFFULL
#define N_CU_OPT1       4

struct ht_entry_opt1_t {
    uint64_t key;
    uint32_t count;
    uint32_t pad;
};

extern "C" void kmer_count_opt1(
    const uint64_t*    kmers,     /* CPU pre-encoded k-mer hashes for this CU */
    ht_entry_opt1_t*   ht_t1,    /* hash table T1 (on-device init) */
    ht_entry_opt1_t*   ht_t2,    /* hash table T2 (on-device init) */
    uint64_t           n_kmers,  /* number of k-mers for this CU */
    uint64_t*          n_unique  /* output */
);
