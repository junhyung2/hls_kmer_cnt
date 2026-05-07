#!/bin/bash
# ============================================================
# 실험: 05_fpga_opt1
# 목적: 3세대 FPGA — CPU pre-encoding + on-device HT init (4CU, 1× U50)
# ============================================================

# [빌드 환경] platform/build_env.txt 참조
# Tool    : Vitis 2022.1, g++ 9.4.0, XRT 2.13.466 (local)
# Platform: xilinx_u50_gen3x16_xdma_5_202210_1

# [실행 환경] platform/run_env.txt 참조
# FPGA    : Alveo U50, 0000:12:00.1 (device 0), Shell UUID 44654095-25B4-C06A-EC6D-0B479D3FEBE8

# [설계 특징]
# - CPU pre-encoding: 호스트에서 FASTQ → canonical k-mer hash 추출
#   → 4개 버킷(h & 3)으로 라우팅 → CU별 uint64_t 배열만 PCIe 전송
#   → seq1/seq2 PCIe 전송 제거, FPGA parse 단계 제거
# - On-device HT init: 커널 내 init_t1()/init_t2() II=1 루프 → DATAFLOW 병렬 초기화
#   → HT 4GB PCIe 전송 제거 (~12s 절감)
# - Cuckoo hash: opt1 동일 (T1/T2, double collision drop)
# - HT_BITS=26 → T1=512MB, T2=512MB per CU × 4CU = 4GB (HBM 상에서 초기화)
# - 설정파일: cfg/connectivity_hw_opt1.cfg (16 HMSS connections)

# [빌드 명령] (로컬에서 완료)
# cd /path/to/src_kmer_opt1
# make -f Makefile_opt1 TARGET=hw kernel   # ~68분 소요
# make -f Makefile_opt1 TARGET=hw host
# scp kmer_count_opt1_hw.xclbin host_opt1_hw u50-server:~/kmer_opt1/

# [실행 명령]
source /opt/xilinx/xrt/setup.sh
cd ~/kmer_opt1

{ time ./host_opt1_hw kmer_count_opt1_hw.xclbin \
    ~/kmer_count/750000read_1.fastq \
    ~/kmer_count/750000read_2.fastq \
    21 ; } 2>&1

# [출력 해석] — host 내 chrono::high_resolution_clock 측정값
# "CPU extraction + file load : 0.83s"
#     t_load0(파일 읽기 시작) → t_load1(extract_kmers_fastq 완료)
# "k-mer PCIe sync : 0.59s"
#     t_sync0 → t_sync1: bo_km[0..3].sync(TO_DEVICE) + bo_nu[0..3].sync(TO_DEVICE)
# "Kernel (HT init + insert) : 3.79s"
#     t0(krnl() 호출) → t1(run.wait() 완료)
# "Total end-to-end : 15.09s"
#     t1 - t_load0  (암묵적 포함: xclbin 로드 + HBM BO 할당 ~9.9s)
