#pragma once
#include <stdint.h>
#ifdef __VITIS_HLS__
#include <ap_int.h>
#endif

#define FE_K                21
#define FE_N_CU             8
#define FE_N_CU_BITS        3
#define FE_LOG2_MAX_PER_CU  23      /* 8M k-mers per CU = 64MB/CU */
#define FE_MAX_PER_CU       (1u << FE_LOG2_MAX_PER_CU)

#define PS_TYPE_PACKED      0x8u
#define PS_TYPE_BOUNDARY    0xBu
#define PS_TYPE_TOKEN       0xFu
#define FE_TOK_RESET        0xFFu
#define FE_TOK_EOF          0xFEu
