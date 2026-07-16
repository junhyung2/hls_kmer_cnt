#!/usr/bin/env bash
# exp73 — build + run the final host. Recovered 2026-07-14: the exact build
# command had never been recorded anywhere, and reconstructing it cost time
# (-luuid is required, and XRT's setup.sh must be sourced or the binary dies
# with "libxrt_coreutil.so.2: cannot open shared object file").
#
# Run this ON THE SERVER (u50-server), from kmer_exp69_h1h2/src.
set -euo pipefail

SRC=exp73_slim.cpp
BIN=h73_slim

XCL=/home/peter04771/kmer_km_single_v15/2xor_11bit_hw.xclbin   # same xclbin as exp50
FQ1=/home/peter04771/kmer_count/ERR022075_1.fastq
FQ2=/home/peter04771/kmer_count/ERR022075_2.fastq

source /opt/xilinx/xrt/setup.sh

g++ -O3 -march=native -std=c++17 \
    -DUSE_REHASH=1 -DBIN_BITS_CFG=12 \
    -I/opt/xilinx/xrt/include \
    "$SRC" -o "$BIN" \
    -L/opt/xilinx/xrt/lib -lxrt_coreutil -luuid -lpthread

# Champion config. -n 3 reports the median e2e (this is the number we quote).
# WALL is reported separately; note -n 3 does three full runs, so for a wall
# figure comparable to KMC you must use -n 1.
/usr/bin/time -f "WALL %e s" ./"$BIN" "$XCL" "$FQ1" "$FQ2" 21 -b 40 -n 3 -as 4

# Expected (ERR022075, 2026-07-10 / re-confirmed 2026-07-14):
#   Unique (acc)              : 109533748   <- must equal KMC exactly
#   dropped inserts (bin full): 0           <- must be 0
#   full bins / total bins    : 0 / 32768
#   occupied slots / capacity : 109533748 / 134217728  (load 81.609%)
#   Median e2e over 3 runs    : ~10.22 s
