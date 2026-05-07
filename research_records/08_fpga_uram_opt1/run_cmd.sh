#!/bin/bash
# ============================================================
# 실험: 08_fpga_uram_opt1
# 목적: 6세대 FPGA — Lazy-clear URAM + U50 × 3카드 12CU 병렬
# ============================================================

# [설계 동기 및 핵심 개선점]
#
# 07_fpga_uram의 병목 분석:
#   kernel 0.290s 중 분해 (300MHz 기준):
#     clear  : 2×131072 cycles × 256 seg = 67M cycles ≈ 224ms  ← 77% 차지
#     insert : 30M k-mer/CU × II=1       = 30M cycles ≈  99ms
#
# 개선 1: Lazy-clear (per-segment clear 제거)
#   각 ht_entry에 uint16_t seg_id 추가.
#   슬롯의 seg_id ≠ 현재 seg → stale(empty로 취급), 덮어씀.
#   세그먼트 루프 전 1회 초기화(262K cycles)만 수행.
#   → clear 오버헤드 67M → 0.26M cycles (257× 감소)
#
# 개선 2: 3카드 × 4CU = 12CU 병렬
#   카드 라우팅: card_id = (h >> 27) % 3
#     비트 [26:10] = T1 인덱스 (커널 내부 사용)
#     비트 [41:27] = 카드 라우팅 → T1 인덱스와 비트 겹침 없음
#   PCIe sync: 카드별 병렬 스레드로 수행
#   각 카드의 kernel이 동시에 실행 (1/3 k-mer씩)
#
# [성능 예측]
#   lazy-clear 후 per-CU cycles: 262K(init) + 30M/3(insert) = 10.26M
#   → kernel ≈ 10.26M / 300MHz ≈ 34ms  (실측 54ms)
#   실측과 예측 차이: URAM init 2-cycle latency, AXI burst 오버헤드 등

# [빌드 환경]
# Tool    : Vitis 2022.1, g++ (XRT 2022.1), XRT 2.13.466
# Platform: xilinx_u50_gen3x16_xdma_5_202210_1
# 빌드 시간: ~59분 (HLS 17분 + impl 42분)
# WNS     : -1.695ns (타이밍 미달이나 실측 정상 동작 확인)

# [실행 환경]
# FPGA: Alveo U50 × 3
#   dev 0: 0000:01:00.1
#   dev 1: 0000:12:00.1
#   dev 2: 0000:17:00.1

# [빌드 명령]
# cd /path/to/src_kmer_sram_opt1
# make -f Makefile_sram_opt1 TARGET=hw kernel   # ~60분
# make -f Makefile_sram_opt1 TARGET=hw host
# scp kmer_count_sram_opt1_hw.xclbin host_sram_opt1_hw u50-server:~/kmer_sram_opt1/

# [실행 명령]
source /opt/xilinx/xrt/setup.sh
mkdir -p ~/kmer_sram_opt1
cd ~/kmer_sram_opt1

{ time ./host_sram_opt1_hw kmer_count_sram_opt1_hw.xclbin \
    ~/kmer_count/750000read_1.fastq \
    ~/kmer_count/750000read_2.fastq \
    21 -d0 0 -d1 1 -d2 2 ; } 2>&1

# [실측 타이밍] (data/output_750k_k21.txt 참조)
# CPU extraction + file load      : 1.072s
# k-mer PCIe sync (parallel)      : 0.355s
# Kernel (3 cards × 4 CU)         : 0.054s   ← 07 대비 5.4× 개선
# Total end-to-end (chrono)       : 20.5s    ← cold XRT 3 xclbin 로드
# Throughput                      : 17,992 MB/s
