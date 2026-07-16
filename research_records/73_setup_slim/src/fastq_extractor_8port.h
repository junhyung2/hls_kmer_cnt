#pragma once
#include <stdint.h>
/* ap_int.h only needed in HLS kernel, not in host g++ compilation */
#ifdef __VITIS_HLS__
#include <ap_int.h>
#endif

/* exp 37: 8 independent AXI ports (one per CU) for km_base
 *
 * Root cause of exp 36 KA=2.028s (expected 0.267s):
 *   Single AXI port → 8 CU regions alternating → no burst coalescing
 *   → single-beat write throughput ≈ 27M writes/s → II_eff ≈ 11-15
 *
 * Fix: km_base_cu0..cu7 each get a separate AXI master (bundle=gmem_km{i})
 *   Each port receives only 1/8 of writes at loop rate 150M writes/s
 *   → per-port write rate = 18.75M/s < per-port AXI throughput (3 outstanding × 7.5M acks/s = 22.5M/s)
 *   → no stall → effective II=2 → KA=0.267s (target)
 *
 * Additionally: host uses parallel radix sort instead of std::sort
 *   std::sort: 960MB working set >> cache → 1.012s
 *   radix sort: O(n) passes, ~0.15-0.2s expected
 */

#define FE_K                21
#define FE_N_CU             8
#define FE_N_CU_BITS        3
#define FE_LOG2_MAX_PER_CU  23      /* 8M k-mers per CU = 64MB per CU (matches 2xor_11bit_hw.xclbin) */
#define FE_MAX_PER_CU       (1u << FE_LOG2_MAX_PER_CU)

/* packed_stream word encoding (upper 4 bits of ap_uint<72>) */
#define PS_TYPE_PACKED      0x8u
#define PS_TYPE_BOUNDARY    0xBu
#define PS_TYPE_TOKEN       0xFu

/* Token sentinels in tok_stream */
#define FE_TOK_RESET        0xFFu
#define FE_TOK_EOF          0xFEu
