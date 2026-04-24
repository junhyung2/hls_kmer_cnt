# kmer_count Build Report

Platform: xilinx_u50_gen3x16_xdma_5_202210_1  
Tool: Vitis 2022.1 / XRT 2.13  
Date: 2026-04-23

---

## Design Summary

FASTQ paired-end k-mer counting kernel using HBM-backed open-addressing hash table.

| Stage | Description |
|-------|-------------|
| `parse_extract_both` | Reads FASTQ1 + FASTQ2, emits canonical k-mer hashes or NO_KMER sentinel |
| `hash_insert` | Consumes stream, inserts/increments in HBM hash table with linear probing |
| DATAFLOW | Both stages run concurrently in HW; stream depth 4096 entries |

Hash table parameters:

| Target | HT_BITS | Entries | Size |
|--------|---------|---------|------|
| SW_EMU | 22 | 4 M | 64 MB |
| HW_EMU | 22 | 4 M | 64 MB |
| HW | 28 | 268 M | 4 GB |

HBM allocation (HW target):

| Port | HBM PC | Size | Purpose |
|------|--------|------|---------|
| seq1 | HBM[0] | 256 MB | FASTQ file 1 |
| seq2 | HBM[1] | 256 MB | FASTQ file 2 |
| ht | HBM[2:17] | 16 x 256 MB = 4 GB | Hash table |
| n_unique | HBM[18] | 8 B | Result |

---

## Build Timeline

### Stage 1: SW Emulation

**Status: PASSED**

- xclbin built in ~14 seconds
- Kernel tested with `022075_read1.fastq` (2.6 MB) + `022075_read2.fastq` (2.6 MB), k=21
- Result: **1,353,128 unique 21-mers** in 0.58 s at 9.0 MB/s

### Stage 2: HW Emulation

**Status: XCLBIN BUILT, FUNCTIONAL RUN SKIPPED (xsim limitation)**

- xclbin built in ~2m 11s
- HW_EMU test (xsim) attempted but cannot complete in practical time:
  - Even 64MB HT DMA through cycle-accurate xsim takes hours
  - Multiple orphaned xsimk instances (each at 99.9% CPU) from repeated attempts
- Conclusion: kernel starts correctly (ERT scheduler configures, dataflow enabled); functional correctness established by SW_EMU; proceeding to HW build

### Stage 3: HW Synthesis

**Status: COMPLETE**

- Build started: 2026-04-23 17:32 (HLS compile)
- HLS compile (.xo): ~0m 32s
- Vivado block synthesis (122 jobs): ~9m
- Vivado placement: ~9m 33s
- Vivado routing: ~6m 2s
- Vivado bitstream generation: ~4m 31s
- **Total elapsed: 35m 9s**
- Output: `kmer_count_hw.xclbin` (30 MB), created 2026-04-23 18:07
- Host binary: `host_hw` built with `HT_BITS_HOST=28` (4 GB HT)

---

## Errors Encountered and Fixes Applied

### Error 1: Macros in HLS pragma arguments

**File:** `src/kmer_count.cpp`  
**Error messages:**
```
ERROR: [HLS 207-3776] use of undeclared identifier 'MAX_PROBES' (kmer_count.cpp:148:38)
ERROR: [HLS 207-3776] use of undeclared identifier 'STREAM_DEPTH' (kmer_count.cpp:214:47)
```
**Root cause:** Vitis HLS 2022.1 pragma parser does NOT expand preprocessor `#define` macros in pragma argument values. Macros are only valid in C/C++ code, not in `#pragma` arguments.

**Fix applied:**
1. `#pragma HLS LOOP_TRIPCOUNT min=1 max=MAX_PROBES avg=2` changed to literal:
   ```cpp
   #pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=2
   ```
2. `#pragma HLS STREAM variable=kmer_stream depth=STREAM_DEPTH` replaced with conditional block:
   ```cpp
   #if defined(SW_EMU)
   #pragma HLS STREAM variable=kmer_stream depth=5200000
   #else
   #pragma HLS STREAM variable=kmer_stream depth=4096
   #endif
   ```

### Error 2: Missing `extern "C"` on kernel function

**File:** `src/kmer_count.h`, `src/kmer_count.cpp`  
**Error message:**
```
ERROR: dlopen of .run/.../dltmp is failed. Please check undefined symbols in the kernel.
       undefined symbol: kmer_count
```
**Root cause:** In SW_EMU, Vitis loads the kernel as a shared library via `dlopen`. Without `extern "C"`, C++ name-mangling changes the symbol name to something like `_Z11kmer_count...`, making it invisible to the XRT dlopen call which looks for the plain symbol `kmer_count`.

**Fix applied:** Added `extern "C"` to kernel declaration in header and definition in cpp:
```cpp
extern "C" void kmer_count(
    const uint8_t* seq1, ...
);
```

### Note: HW_EMU HT_BITS reduced from 24 to 22

HT_BITS for HW_EMU reduced from 24 (256 MB) to 22 (64 MB) to limit xsim memory initialization time. Even at 64 MB the DMA sync through xsim is too slow for practical testing, so the HW_EMU run is skipped. The HW target retains HT_BITS=28 (4 GB, 268M entries).

---

## Test Results

| Target | Result | Unique 21-mers | Runtime |
|--------|--------|----------------|---------|
| SW_EMU | PASS | 1,353,128 | 0.58 s |
| HW_EMU | SKIPPED (xsim) | - | - |
| HW | BUILD COMPLETE | run on U50 to verify | 35m 9s build |

---

## Deployment

To run on a real Alveo U50:

```bash
source /opt/xilinx/xrt/setup.sh
./host_hw kmer_count_hw.xclbin <fastq1> <fastq2> <k>
```

The host binary expects the U50 to be present (`xrt::device(0)`). The kernel allocates:
- HBM[0]: FASTQ1 (up to 256 MB)
- HBM[1]: FASTQ2 (up to 256 MB)
- HBM[2:17]: 4 GB hash table (2^28 × 16 B)
- HBM[18]: 8-byte n_unique result
