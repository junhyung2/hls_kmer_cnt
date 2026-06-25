#!/bin/bash
# exp58: sub-bin scatter + 8-CU parallel radix sort
# Run on FPGA server (u50-server)
# xclbin: ~/kmer_sort58/sort58_hw.xclbin
# host:   ~/kmer_sort58/host_sort58

source /opt/xilinx/xrt/setup.sh
cd ~/kmer_sort58

# 21G
./host_sort58 sort58_hw.xclbin \
    ~/kmer_count_err3239334/ERR3239334_21G_1.fastq \
    ~/kmer_count_err3239334/ERR3239334_21G_2.fastq \
    21 -b 64

# 50G
./host_sort58 sort58_hw.xclbin \
    ~/kmer_count_err3239334/ERR3239334_50G_1.fastq \
    ~/kmer_count_err3239334/ERR3239334_50G_2.fastq \
    21 -b 64
