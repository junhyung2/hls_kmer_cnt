#!/bin/bash
# ============================================================
# 실험: 06_fpga_card2
# 목적: 4세대 FPGA — Dual U50 (8CU 병렬, 동일 xclbin × 2카드)
# ============================================================

# [빌드 환경] platform/build_env.txt 참조
# Tool    : Vitis 2022.1, g++ 9.4.0, XRT 2.13.466 (local)
# Platform: xilinx_u50_gen3x16_xdma_5_202210_1
# 주의     : xclbin은 05_fpga_opt1과 동일 파일 — 호스트만 다름

# [실행 환경] platform/run_env.txt 참조
# FPGA    : Alveo U50 × 2
#   Card 0 : 0000:12:00.1, hbm_aclk 450MHz
#   Card 1 : 0000:17:00.1, hbm_aclk 415MHz
#   두 카드 Shell UUID 동일: 44654095-25B4-C06A-EC6D-0B479D3FEBE8

# [설계 특징]
# - 커널: kmer_count_opt1 과 완전 동일 (4CU 설계)
# - 호스트 변경점:
#     · 8버킷 라우팅 (h & 7) → card0=버킷[0-3], card1=버킷[4-7]
#     · xrt::device 2개 독립 open → 각각 xclbin 로드
#     · 4CU × 2카드 = 8 run 동시 launch, 순차 wait
#     · PCIe sync: 8 BO 순차 전송 (병렬화 미적용 → 1.15s)
# - xclbin 파일: kmer_count_opt1_hw.xclbin (opt1과 동일)
# - 설정파일: cfg/connectivity_hw_opt1.cfg (각 카드에 독립 적용)

# [빌드 명령]
# cd /path/to/src_kmer_opt_card2
# make -f Makefile_card2 TARGET=hw host   # xclbin은 opt1 것 재사용
# scp host_card2_hw u50-server:~/kmer_opt1/

# [실행 명령]
source /opt/xilinx/xrt/setup.sh
cd ~/kmer_opt1

{ time ./host_card2_hw kmer_count_opt1_hw.xclbin \
    ~/kmer_count/750000read_1.fastq \
    ~/kmer_count/750000read_2.fastq \
    21 -d0 0 -d1 1 ; } 2>&1

# 옵션: -d0 <device_index_card0>  -d1 <device_index_card1>
# device index 확인: xbutil examine → BDF 순서대로 0, 1

# [출력 해석] — host 내 chrono::high_resolution_clock 측정값
# "CPU extraction + file load : 0.80s"
#     8버킷 라우팅 포함 (h & 7), opt1과 실질적으로 동일
# "k-mer sync (both cards) : 1.15s"
#     t_sync0 → t_sync1: 8 bo_km.sync(TO_DEVICE) 순차 실행
#     (병렬화 시 ~0.6s 예상: std::thread로 각 카드 sync 분리)
# "Kernel (both cards) : 1.96s"
#     t0(8× launch) → t1(8× wait 완료) — 두 카드 물리적 병렬 실행
# "Total end-to-end : 13.17s"
#     t1 - t_load0  (암묵적: 2× xclbin 로드 + 2× HBM BO 할당 ~10.2s)
