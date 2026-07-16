#!/usr/bin/env bash
# exp73 — regenerate the FPGA bitstream (2xor_11bit_hw.xclbin).
#
# exp73 is a HOST-ONLY result: it runs on the SAME xclbin as exp50, built here
# in exp40 (the "2-XOR hash + II=1" kernel). No resynthesis was done from exp50
# through exp73. You only need this if the xclbin is lost — the prebuilt binary
# already lives on the server at:
#     ~/kmer_km_single_v15/2xor_11bit_hw.xclbin
# (that is the exact path build_and_run.sh's XCL points at).
#
# Kernel source md5 (must match, identical across exp40/exp50):
#     c768aaffb4999d1752271fbcc2c05d7f  fastq_extractor_8port_2xor.cpp
#
# Synthesis takes ~2 h. Run ON THE SERVER (u50-server). ~28.7 s of the reported
# "cold" wall is one-time bitstream programming of the 3 cards, not this build.
set -euo pipefail

PLATFORM=xilinx_u50_gen3x16_xdma_5_202210_1

source /tools/Xilinx/Vitis/2023.1/settings64.sh
source /opt/xilinx/xrt/setup.sh

# Step 1: compile the kernel (2-XOR hash, II=1 @ 300 MHz)
v++ -c --platform "$PLATFORM" --target hw \
    -I . -k fastq_extractor_8port \
    -o fastq_extractor_8port_hw.xo \
    fastq_extractor_8port_2xor.cpp

# Step 2: link 3 CUs, 8 km_base HBM banks per CU (see connectivity_hw.cfg)
v++ -l --platform "$PLATFORM" --target hw \
    --config connectivity_hw.cfg \
    -o 2xor_11bit_hw.xclbin \
    fastq_extractor_8port_hw.xo

echo "built 2xor_11bit_hw.xclbin — install it where build_and_run.sh's XCL points."
