/* exp58: fastq_extractor_8port_nohash.cpp
 *
 * vs exp50/55 (2-XOR hash):
 *   exp50/55: output = hash64_fe(canonical)  → hash used only for hash-table routing (exp50)
 *             radix sort (exp55) doesn't need hashing; hash was unnecessary overhead
 *
 *   exp58: output = canonical k-mer (42-bit) directly
 *          CU routing: (y ^ (y >> 21)) & 7  — single XOR-shift, same distribution quality
 *          EXTRACT_LOOP critical path: reduced ~2.3ns → ~1.8ns (more timing margin at 300MHz)
 *
 * Host sorts canonical values directly without hash computation.
 * Interface is identical to exp50/55 (same argument list, same AXI pragmas).
 */

#include "fastq_extractor_8port.h"
#include <hls_stream.h>
#include <ap_int.h>

static inline uint64_t swar_has_newline(uint64_t v)
{
#pragma HLS INLINE
    uint64_t x = v ^ 0x0a0a0a0a0a0a0a0aULL;
    return (x - 0x0101010101010101ULL) & ~x & 0x8080808080808080ULL;
}

static inline uint8_t encode_base(uint8_t c)
{
#pragma HLS INLINE
    if (c == 'A' || c == 'a') return 0;
    if (c == 'C' || c == 'c') return 1;
    if (c == 'G' || c == 'g') return 2;
    if (c == 'T' || c == 't') return 3;
    return 4;
}

static inline int count_newlines(uint64_t nl_mask)
{
#pragma HLS INLINE
    int n = 0;
    for (int b = 0; b < 8; b++) {
#pragma HLS UNROLL
        n += (int)((nl_mask >> (b * 8 + 7)) & 1ULL);
    }
    return n;
}

static void scan_fastq(
    const uint64_t*           fastq_data,
    uint64_t                  fastq_len,
    hls::stream<ap_uint<72>>& packed_out
) {
    int line_num = 0;
    const uint64_t n_chunks = (fastq_len + 7) / 8;

SCAN_LOOP:
    for (uint64_t ci = 0; ci < n_chunks; ci++) {
#pragma HLS PIPELINE II=2
#pragma HLS LOOP_TRIPCOUNT min=4000000 max=33000000 avg=16250000

        uint64_t chunk   = fastq_data[ci];
        uint64_t nl_mask = swar_has_newline(chunk);
        bool     has_nl  = (nl_mask != 0);

        if (!has_nl && line_num != 1) continue;

        ap_uint<72> out_word;
        if (!has_nl) {
            ap_uint<32> pb = 0;
            for (int b = 0; b < 8; b++) {
#pragma HLS UNROLL
                pb.range(b*4+3, b*4) = encode_base((uint8_t)(chunk >> (b*8)));
            }
            out_word = ((ap_uint<72>)PS_TYPE_PACKED << 68) | (ap_uint<72>)pb;
        } else {
            out_word = ((ap_uint<72>)PS_TYPE_BOUNDARY << 68) |
                       ((ap_uint<72>)(line_num & 3) << 64) |
                       (ap_uint<72>)chunk;
            int nl_cnt = count_newlines(nl_mask);
            line_num = (line_num + nl_cnt) & 3;
        }
        packed_out.write(out_word);
    }

    packed_out.write(((ap_uint<72>)PS_TYPE_TOKEN << 68) | (ap_uint<72>)FE_TOK_RESET);
    packed_out.write(((ap_uint<72>)PS_TYPE_TOKEN << 68) | (ap_uint<72>)FE_TOK_EOF);
}

static void unpack_tokens(
    hls::stream<ap_uint<72>>& packed_in,
    hls::stream<uint8_t>&     tok_out
) {
    const uint32_t n_iter = 400000000U;

    ap_uint<72> cur_word = 0;
    ap_uint<4>  cur_type = PS_TYPE_TOKEN;
    uint8_t     byte_pos = 8;
    int         cur_ln   = 0;

    ap_uint<4>  bases[8];
#pragma HLS ARRAY_PARTITION variable=bases complete dim=1
    ap_uint<8>  rbytes[8];
#pragma HLS ARRAY_PARTITION variable=rbytes complete dim=1

UNPACK_FSM:
    for (uint32_t wi = 0; wi < n_iter; wi++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=2000000 max=133000000 avg=63000000

        if (byte_pos >= 8) {
            cur_word = packed_in.read();
            cur_type = cur_word.range(71, 68);
            byte_pos = 0;

            if (cur_type == PS_TYPE_PACKED) {
                for (int b = 0; b < 8; b++) {
#pragma HLS UNROLL
                    bases[b] = cur_word.range(4*b+3, 4*b);
                }
            } else if (cur_type == PS_TYPE_BOUNDARY) {
                cur_ln = (int)(uint8_t)cur_word.range(65, 64);
                for (int b = 0; b < 8; b++) {
#pragma HLS UNROLL
                    rbytes[b] = (ap_uint<8>)cur_word.range(8*b+7, 8*b);
                }
            }
        }

        if (cur_type == PS_TYPE_PACKED) {
            tok_out.write((uint8_t)bases[byte_pos]);
            byte_pos = (byte_pos == 7) ? 8 : (uint8_t)(byte_pos + 1);

        } else if (cur_type == PS_TYPE_BOUNDARY) {
            uint8_t c = (uint8_t)rbytes[byte_pos];
            if (c == '\n') {
                if (cur_ln == 1) tok_out.write(FE_TOK_RESET);
                cur_ln = (cur_ln < 3) ? cur_ln + 1 : 0;
            } else if (cur_ln == 1) {
                tok_out.write(encode_base(c));
            }
            byte_pos = (byte_pos == 7) ? 8 : (uint8_t)(byte_pos + 1);

        } else {
            uint8_t t = (uint8_t)cur_word.range(7, 0);
            tok_out.write(t);
            byte_pos = 8;
            if (t == FE_TOK_EOF) break;
        }
    }
}

/* exp58: no hash — output canonical k-mer (42-bit) directly.
 * CU routing: (y ^ (y >> 21)) & 7  (single XOR-shift, ~0.5ns critical path)
 * vs exp50/55: hash64_fe used two XOR-shifts (~1.0ns) that were unnecessary for radix sort.
 */
static void extract_kmers_8port(
    hls::stream<uint8_t>& seq_in,
    uint64_t* km_base_cu0, uint64_t* km_base_cu1,
    uint64_t* km_base_cu2, uint64_t* km_base_cu3,
    uint64_t* km_base_cu4, uint64_t* km_base_cu5,
    uint64_t* km_base_cu6, uint64_t* km_base_cu7,
    uint64_t* se_out
) {
    uint64_t flat_cnt[FE_N_CU] = {};
#pragma HLS ARRAY_PARTITION variable=flat_cnt complete dim=1

    const uint64_t km_mask = (1ULL << (2 * FE_K)) - 1;
    const int      shift   = (FE_K - 1) * 2;

    uint64_t kmer_fwd = 0;
    uint64_t kmer_rev = 0;
    int      kmer_len = 0;

    const uint32_t n_tokens = 500000000U;

EXTRACT_LOOP:
    for (uint32_t ti = 0; ti < n_tokens; ti++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=333000 max=166000000 avg=65000000
#pragma HLS DEPENDENCE variable=flat_cnt inter false

        uint8_t tok = seq_in.read();

        if (tok == FE_TOK_EOF)   break;
        if (tok == FE_TOK_RESET) { kmer_fwd = 0; kmer_rev = 0; kmer_len = 0; continue; }
        if (tok == 4)            { kmer_fwd = 0; kmer_rev = 0; kmer_len = 0; continue; }

        kmer_fwd = ((kmer_fwd << 2) | (uint64_t)tok) & km_mask;
        kmer_rev = (kmer_rev >> 2) | ((uint64_t)(3 - tok) << shift);
        kmer_len++;

        if (kmer_len >= FE_K) {
            uint64_t y      = (kmer_fwd < kmer_rev) ? kmer_fwd : kmer_rev;
            int      cu     = (int)((y ^ (y >> 21)) & (uint64_t)(FE_N_CU - 1));
            uint64_t c      = flat_cnt[cu];

            switch (cu) {
            case 0: km_base_cu0[c] = y; break;
            case 1: km_base_cu1[c] = y; break;
            case 2: km_base_cu2[c] = y; break;
            case 3: km_base_cu3[c] = y; break;
            case 4: km_base_cu4[c] = y; break;
            case 5: km_base_cu5[c] = y; break;
            case 6: km_base_cu6[c] = y; break;
            case 7: km_base_cu7[c] = y; break;
            }
            flat_cnt[cu] = c + 1;
        }
    }

    for (int cu = 0; cu < FE_N_CU; cu++) {
#pragma HLS PIPELINE II=1
        se_out[cu] = flat_cnt[cu];
    }
}

extern "C" void fastq_extractor_8port(
    const uint64_t* fastq_data,
    uint64_t        fastq_len,
    uint64_t*       km_base_cu0,
    uint64_t*       km_base_cu1,
    uint64_t*       km_base_cu2,
    uint64_t*       km_base_cu3,
    uint64_t*       km_base_cu4,
    uint64_t*       km_base_cu5,
    uint64_t*       km_base_cu6,
    uint64_t*       km_base_cu7,
    uint64_t*       se_out
) {
#pragma HLS INTERFACE m_axi port=fastq_data bundle=gmem_fq \
        max_read_burst_length=16 num_read_outstanding=8 \
        max_write_burst_length=2  num_write_outstanding=1

#pragma HLS INTERFACE m_axi port=km_base_cu0 bundle=gmem_km0 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu1 bundle=gmem_km1 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu2 bundle=gmem_km2 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu3 bundle=gmem_km3 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu4 bundle=gmem_km4 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu5 bundle=gmem_km5 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu6 bundle=gmem_km6 \
        max_write_burst_length=16 num_write_outstanding=128
#pragma HLS INTERFACE m_axi port=km_base_cu7 bundle=gmem_km7 \
        max_write_burst_length=16 num_write_outstanding=128

#pragma HLS INTERFACE m_axi port=se_out bundle=gmem_se \
        max_write_burst_length=8 num_write_outstanding=1

#pragma HLS INTERFACE s_axilite port=fastq_data   bundle=control
#pragma HLS INTERFACE s_axilite port=fastq_len    bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu0  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu1  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu2  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu3  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu4  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu5  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu6  bundle=control
#pragma HLS INTERFACE s_axilite port=km_base_cu7  bundle=control
#pragma HLS INTERFACE s_axilite port=se_out        bundle=control
#pragma HLS INTERFACE s_axilite port=return        bundle=control

    hls::stream<ap_uint<72>> packed_stream("packed_stream");
#pragma HLS STREAM variable=packed_stream depth=4096
    hls::stream<uint8_t> tok_stream("tok_stream");
#pragma HLS STREAM variable=tok_stream depth=32768

#pragma HLS DATAFLOW
    scan_fastq(fastq_data, fastq_len, packed_stream);
    unpack_tokens(packed_stream, tok_stream);
    extract_kmers_8port(tok_stream,
        km_base_cu0, km_base_cu1, km_base_cu2, km_base_cu3,
        km_base_cu4, km_base_cu5, km_base_cu6, km_base_cu7,
        se_out);
}
