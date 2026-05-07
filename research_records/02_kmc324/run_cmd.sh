#!/bin/bash
# ============================================================
# 실험: 02_kmc324
# 목적: 업계 표준 CPU k-mer 카운팅 도구 (멀티스레드 기준선)
# ============================================================

# [실행 환경]
# 머신     : u50-server (AMD Ryzen 7 9800X3D, 16 logical cores, Ubuntu 20.04.6)
# 상세     : platform/run_env.txt 참조
# 바이너리  : KMC 3.2.4 정적 링크 바이너리 (소스 불포함)
#   다운로드 : https://github.com/refresh-bio/KMC/releases/tag/v3.2.4
#   파일명   : KMC-3.2.4/bin/kmc  (Linux x86_64 정적 바이너리, 추가 의존성 없음)

# [데이터셋]
# 파일     : platform/dataset.txt 참조
# 경로     : ~/kmer_count/750000read_1.fastq, ~/kmer_count/750000read_2.fastq

# [실행 명령]
cd ~/kmer_count
mkdir -p /tmp/kmc_tmp

{ time ~/kmc \
    -k21 \
    -t4 \
    -ci1 \
    @<(echo 750000read_1.fastq; echo 750000read_2.fastq) \
    /tmp/kmc_out \
    /tmp/kmc_tmp ; } 2>&1

rm -rf /tmp/kmc_out* /tmp/kmc_tmp

# [옵션 설명]
# -k21   : k-mer 길이 21
# -t4    : 4 스레드
# -ci1   : count ≥ 1 인 k-mer만 포함 (= 전체 unique)
# @<(...): 입력 파일 목록 (process substitution)
# /tmp/kmc_out    : 출력 파일 prefix
# /tmp/kmc_tmp    : 임시 디렉토리

# [출력 해석]
# stdout: "No. of unique k-mers : 9,044,812"
# stderr (bash time):
#   real = 2.107s (벽시계), user = 4.472s (4스레드 합산 CPU time)
