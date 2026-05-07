#!/bin/bash
# ============================================================
# 실험: 07_fpga_uram
# 목적: 5세대 FPGA — Alveo U50 URAM을 hash table로 활용해
#        HBM 랜덤 접근 latency 병목 제거 (1카드 4CU)
# ============================================================

# [설계 동기]
# kc-c4의 sub hash table 방식(전체 HT를 소규모 sub-table로 분할 →
# 각 sub-table이 L3 cache에 수용 → cache miss 제거)을 FPGA에 적용:
#   • HBM (opt1/card2 기존 설계) : 랜덤 접근 latency ~100ns
#   • URAM (본 설계)             : 접근 latency ~2 cycles ≈ 6ns  (50× 개선)
#
# [핵심 파라미터]
#   N_CU   = 4   (카드당 CU 수, opt1과 동일)
#   N_SEG  = 16  (CU당 URAM 세그먼트 수)
#   SRAM_HT_BITS = 18  → T_SIZE = 131072 entries per T per segment
#   T1+T2 size per CU = 2 × 2MB = 4MB → 128 URAM / CU × 4CU = 512 URAM (80% 사용)
#
# [부하율 검증]
#   unique k-mer / segment / CU = 9M / (4 × 16) ≈ 141K
#   T1+T2 slots                  = 2 × 131072 = 262K
#   load factor                  ≈ 0.54  (OK)
#
# [HBM 사용량 비교]
#   opt1 : kmers(960MB) + T1(2GB) + T2(2GB) = ~5GB BO 할당
#   sram : kmers(960MB) + seg_ends(512B) + n_unique(32B) → ~960MB only
#   → XRT BO 할당 시간 대폭 감소 (~6s → ~1s)

# [빌드 환경] platform/build_env.txt 참조
# Tool    : Vitis 2022.1, g++ 9.4.0, XRT 2.13.466 (local)
# Platform: xilinx_u50_gen3x16_xdma_5_202210_1
# 빌드 디렉토리: src_kmer_sram/

# [실행 환경] platform/run_env.txt 참조
# FPGA    : Alveo U50, 0000:12:00.1 (device 0)

# [빌드 명령]
# cd /path/to/src_kmer_sram
# make -f Makefile_sram TARGET=hw kernel   # ~60~90분 소요
# make -f Makefile_sram TARGET=hw host
# scp kmer_count_sram_hw.xclbin host_sram_hw u50-server:~/kmer_sram/

# [실행 명령]
source /opt/xilinx/xrt/setup.sh
mkdir -p ~/kmer_sram
cd ~/kmer_sram

{ time ./host_sram_hw kmer_count_sram_hw.xclbin \
    ~/kmer_count/750000read_1.fastq \
    ~/kmer_count/750000read_2.fastq \
    21 -d 0 ; } 2>&1

# [예상 타이밍] (이론값, 300MHz 기준)
# URAM clear    : 2 × 131072 cycles × 16 seg / 300M ≈ 14 ms/CU
# URAM insert   : ~30M k-mer × II≈1 / 300M    ≈ 100 ms/CU
# Kernel (4CU) : ≈ 114 ms   (vs opt1 3790 ms = 33× speedup)
# XRT overhead  : ~4s        (xclbin 로드 + 960MB BO 할당; opt1 9.9s 대비 ~2.5× 감소)
# Total e2e     : ~5.5s      (vs opt1 15.1s = 2.7× speedup)
#
# 실측값은 data/output_750k_k21.txt 참조
