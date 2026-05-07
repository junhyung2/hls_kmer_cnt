#!/bin/bash
# ============================================================
# 실험: 03_fpga_original
# 목적: 1세대 FPGA 설계 — 단일 AXI master, linear probe loop
# ============================================================

# [빌드 환경] platform/build_env.txt 참조
# Tool    : Vitis 2022.1, g++ 9.4.0, XRT 2.13.466 (local)
# Platform: xilinx_u50_gen3x16_xdma_5_202210_1

# [실행 환경] platform/run_env.txt 참조
# 머신    : u50-server (AMD Ryzen 7 9800X3D, RAM 186GB)
# XRT     : 2.15.225 (2023.1)
# FPGA    : Alveo U50, 0000:17:00.1 (device 1), Shell UUID 44654095-25B4-C06A-EC6D-0B479D3FEBE8

# [설계 특징]
# - HT_BITS=26 → T_SIZE=2^25=32M entries per T, 1 table, HT 1GB/CU × 4CU = 4GB
# - 단일 AXI master (gmem2): T1/T2 구분 없이 하나의 HT
# - Collision 해결: linear probe, 최대 128 스텝 → 외부 루프 II>>1
# - 호스트가 4GB HT를 memset(EMPTY) 후 PCIe로 전송 → DMA 비용 큼
# - 설정파일: cfg/connectivity_hw_opt.cfg

# [빌드 명령] (서버에서 이미 완료, 재빌드 시 참고)
# cd ~/kmer_count
# v++ -c --platform xilinx_u50_gen3x16_xdma_5_202210_1 --target hw \
#     -I src -k kmer_count_opt -o kmer_count_opt_hw.xo src/kmer_count_opt.cpp
# v++ -l --platform xilinx_u50_gen3x16_xdma_5_202210_1 --target hw \
#     --config connectivity_hw_opt.cfg -o kmer_count_opt_hw.xclbin kmer_count_opt_hw.xo
# g++ -std=c++17 -O2 -DHT_BITS_HOST=26 -I/opt/xilinx/xrt/include -I src \
#     -o host_opt_hw src/host_opt.cpp -L/opt/xilinx/xrt/lib -lxrt_coreutil -lpthread

# [실행 명령]
source /opt/xilinx/xrt/setup.sh
cd ~/kmer_count

./host_opt_hw kmer_count_opt_hw.xclbin \
    750000read_1.fastq \
    750000read_2.fastq \
    21 -d 1

# [출력 해석]
# stdout : "All CUs finished in 19.7492 s" = 커널 실행 시간
# data/xrt_summary.csv : APPLICATION_RUN_TIME_MS = 30392.2 → 총 30.4s
#   (xclbin 로드 + 4GB HT memset + PCIe DMA + 커널 + 결과 회수 포함)
