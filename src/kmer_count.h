#ifndef KMER_COUNT_H
#define KMER_COUNT_H

#include <stdint.h>
#include <hls_stream.h>
#include <ap_int.h>

/* Hash table sizing per build target
 * SW_EMU  : 2^22 =   4 M entries x 16 B =   64 MB  (1 HBM PC, fast sim)
 * HW_EMU  : 2^24 =  16 M entries x 16 B =  256 MB  (1 HBM PC)
 * HW      : 2^28 = 268 M entries x 16 B = 4096 MB  (16 HBM PCs, 400 MB FASTQ)
 */
#if defined(SW_EMU)
  #define HT_BITS  22
#elif defined(HW_EMU)
  #define HT_BITS  22
#else
  #define HT_BITS  28
#endif

#define HT_SIZE    (1ULL << HT_BITS)
#define HT_MASK    (HT_SIZE - 1ULL)

/* 0xFFFF...FF cannot be hash64 output for k<=31
 * (hash64 output is always masked to 2^(2k)-1 < 2^62) */
#define EMPTY_KEY  0xFFFFFFFFFFFFFFFFULL
/* same sentinel used as "no k-mer produced for this byte" in the stream */
#define NO_KMER    0xFFFFFFFFFFFFFFFFULL

/* Maximum linear-probe steps before giving up.
 * At 50% load factor P(probes > 128) < 2^-128: effectively impossible. */
#define MAX_PROBES 128

/* Stream depth between parse_extract and hash_insert:
 *   SW_EMU  - DATAFLOW runs sequentially in simulation; stream must buffer
 *             the full output of parse_extract (2 x 2.5 MB test files ->
 *             at most ~5 M stream entries).
 *   HW/HW_EMU - true concurrent execution; 4 K entries is sufficient. */
#if defined(SW_EMU)
  #define STREAM_DEPTH  5200000
#else
  #define STREAM_DEPTH  4096
#endif

#define MAX_K  31
#define MIN_K   1

/* Hash table entry: 16 B (key 8 B + count 4 B + pad 4 B).
 * Aligned to AXI 128-bit beat; two entries per 256-bit HBM access. */
typedef struct {
    uint64_t key;    /* EMPTY_KEY = vacant slot */
    uint32_t count;
    uint32_t pad;
} ht_entry_t;

/* Top-level HLS kernel
 * seq1     (m_axi) : FASTQ file 1  -> HBM[0]
 * seq2     (m_axi) : FASTQ file 2  -> HBM[1]
 * ht       (m_axi) : hash table    -> HBM[2]     (SW/HW_EMU)
 *                                  -> HBM[2:17]  (HW, 4 GB)
 * n_unique (m_axi) : unique k-mer count output -> HBM[18] (HW) */
extern "C" void kmer_count(
    const uint8_t* seq1,
    const uint8_t* seq2,
    ht_entry_t*    ht,
    uint64_t       len1,
    uint64_t       len2,
    int            k,
    uint64_t*      n_unique
);

#endif /* KMER_COUNT_H */
