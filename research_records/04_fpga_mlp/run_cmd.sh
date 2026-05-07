#!/bin/bash
# ============================================================
# 실험: 04_fpga_mlp
# 목적: 2세대 FPGA — MLP (Memory-Level Parallelism), cuckoo hash, DATAFLOW
# ============================================================

# [빌드 환경] platform/build_env.txt 참조
# Tool    : Vitis 2022.1, g++ 9.4.0, XRT 2.13.466 (local)
# Platform: xilinx_u50_gen3x16_xdma_5_202210_1

# [실행 환경] platform/run_env.txt 참조
# FPGA    : Alveo U50, 0000:12:00.1 (device 0), Shell UUID 44654095-25B4-C06A-EC6D-0B479D3FEBE8

# [설계 특징]
# - HT_BITS=26 → T_SIZE=2^25=32M entries, T1+T2 = 512MB+512MB = 1GB/CU × 4CU = 4GB
# - AXI master 2개 분리: gmem_t1(T1 전용) / gmem_t2(T2 전용) → 동시 발행 가능
# - Cuckoo hash: 직접 탐색(T1[h1], T2[h2]) → double collision 시 drop
# - PIPELINE II=1 + 32 outstanding requests → HBM ~30cycle 레이턴시 은닉
# - DATAFLOW: parse_extract ↔ cuckoo_insert 병렬 실행 (stream 통신)
# - 호스트: seq1/seq2 PCIe 전송 + HT 4GB memset 후 PCIe 전송 (주 병목)
# - 설정파일: cfg/connectivity_hw_mlp.cfg (20 HMSS connections)

# [빌드 명령] (로컬에서 완료, 재빌드 시)
# cd /path/to/src_kmer_mlp
# make -f Makefile_mlp TARGET=hw kernel   # ~70분 소요

# [실행 명령]
source /opt/xilinx/xrt/setup.sh
cd ~/kmer_mlp

{ time ./host_mlp kmer_mlp.xclbin \
    ~/kmer_count/750000read_1.fastq \
    ~/kmer_count/750000read_2.fastq \
    21 ; } 2>&1

# [출력 해석]
# "[kmer_count_mlp] All CUs finished in 5.25093 s" = 커널 시간
# bash time real = 17.98s = xclbin 로드 + HT 4GB DMA(~12.4s) + 커널(5.25s) + 기타
