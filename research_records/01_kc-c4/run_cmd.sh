#!/bin/bash
# ============================================================
# 실험: 01_kc-c4
# 목적: CPU 기준 k-mer 카운팅 (단일스레드 기준선)
# ============================================================

# [실행 환경]
# 머신     : u50-server (AMD Ryzen 7 9800X3D, RAM 186GB, Ubuntu 20.04.6)
# 상세     : platform/run_env.txt 참조
# 바이너리  : src/kc-c4.c 를 서버에서 직접 컴파일
#   빌드 명령: gcc -O2 -o kc-c4 kc-c4.c -lpthread
#   (원본 kseq.h 포함 필요; zlib 의존성 제거 패치 적용됨 — src/kc-c4.c 참조)

# [데이터셋]
# 파일     : platform/dataset.txt 참조
# 경로     : ~/kmer_count/750000read_1.fastq, ~/kmer_count/750000read_2.fastq
# MD5      : dc419fd9... / f061015f...

# [실행 명령]
cd ~/kmer_count

# kc-c4는 단일 파일만 입력받으므로 두 파일 합산 후 실행
cat 750000read_1.fastq 750000read_2.fastq > /tmp/combined_750k.fastq
{ time ~/kmer-cnt/kc-c4 /tmp/combined_750k.fastq 21 ; } 2>&1
rm -f /tmp/combined_750k.fastq

# [출력 해석]
# stdout: "count<TAB>frequency" 형식 히스토그램
#   → 1열 합산 = total unique k-mers = 9,044,812
# stderr (bash time):
#   real  = 벽시계 기준 총 경과 시간 → 0.578s
#   user  = CPU 사용 시간 합산 (멀티코어 효과 반영) → 2.323s
#   sys   = 시스템콜 시간
