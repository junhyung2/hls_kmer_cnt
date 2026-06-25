#!/bin/bash
set -e
cd ~/kmer_sort58
source /tools/Xilinx/Vitis/2023.1/settings64.sh
source /opt/xilinx/xrt/setup.sh

PLATFORM=xilinx_u50_gen3x16_xdma_5_202210_1
LOG=build_hw.log

echo "[$(date)] BUILD START (exp58: no-hash kernel + sub-bin scatter sort)" | tee $LOG

echo "[$(date)] Step 1: v++ compile (HLS → XO)" | tee -a $LOG
v++ -c --platform $PLATFORM --target hw \
    -I src \
    -k fastq_extractor_8port \
    -o fastq_extractor_8port_nohash.xo \
    src/fastq_extractor_8port_nohash.cpp 2>&1 | tee -a $LOG
echo "[$(date)] Step 1 DONE" | tee -a $LOG

echo "[$(date)] Step 2: v++ link (XO → xclbin, 300MHz)" | tee -a $LOG
v++ -l --platform $PLATFORM --target hw \
    --config cfg/connectivity_58.cfg \
    --kernel_frequency 300 \
    -o sort58_hw.xclbin \
    fastq_extractor_8port_nohash.xo 2>&1 | tee -a $LOG
echo "[$(date)] Step 2 DONE" | tee -a $LOG

echo "[$(date)] Step 3: compile host" | tee -a $LOG
g++ -O3 -std=c++17 \
    -I src \
    -I /opt/xilinx/xrt/include \
    -L /opt/xilinx/xrt/lib \
    -o host_sort58 \
    src/host_3card_sort58.cpp \
    -lxrt_coreutil -lpthread -luuid 2>&1 | tee -a $LOG
echo "[$(date)] Step 3 DONE" | tee -a $LOG

echo "[$(date)] BUILD COMPLETE" | tee -a $LOG
