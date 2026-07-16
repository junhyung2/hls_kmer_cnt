# exp73 — final result (self-contained reproduction bundle)

exp73 is the champion result of this project. On ERR022075 (22.72M read pairs,
k=21) it beats KMC 3.2.4 `-t16` on **both** e2e (1.39×) and full wall (1.05×),
with **identical** unique count (109,533,748) and **zero resynthesis** from exp50.

Every gain from exp50 → exp73 is host-side. exp73 runs on the *same* bitstream
as exp50 (`2xor_11bit_hw.xclbin`, the exp40 "2-XOR hash + II=1" kernel).

## What's in this bundle

Before this cleanup the record held only the host source, the host build/run
script, and the analysis — you could not rebuild the host (its header was
missing) nor regenerate the bitstream (kernel sources were elsewhere). Now the
record is self-contained:

```
73_setup_slim/
├── README.md                       ← this file
├── build_and_run.sh                host build + champion run (server: kmer_exp69_h1h2/src)
├── src/
│   ├── exp73_slim.cpp              the host (3× U50, host-side counting)
│   └── fastq_extractor_8port.h     shared kernel/host header — REQUIRED to compile the host
├── kernel/                         everything to regenerate the xclbin (only if it's lost)
│   ├── fastq_extractor_8port_2xor.cpp   kernel source (md5 c768aaffb4999d1752271fbcc2c05d7f)
│   ├── fastq_extractor_8port.h          same shared header
│   ├── connectivity_hw.cfg              3 CUs, per-CU HBM bank mapping
│   └── build_xclbin.sh                  v++ compile + link (~2 h synthesis)
└── data/
    └── analysis.txt                the 2026-07-14 comment-correction writeup
```

Not included (binaries, live on the server, too large / machine-specific):
- `2xor_11bit_hw.xclbin` — the prebuilt bitstream, at `~/kmer_km_single_v15/2xor_11bit_hw.xclbin`
- `h73_slim` — the built host binary, at `~/kmer_exp69_h1h2/src/h73_slim`
- ERR022075 FASTQ, at `~/kmer_count/ERR022075_{1,2}.fastq`

## How to reproduce

1. **Host only** (the normal case — the xclbin already exists on the cards):
   run `build_and_run.sh` on u50-server from `~/kmer_exp69_h1h2/src`. It compiles
   with `-DUSE_REHASH=1 -DBIN_BITS_CFG=12` and runs the champion config
   `... 21 -b 40 -n 3 -as 4`. Needs `fastq_extractor_8port.h` next to the source
   (now in `src/`), `-luuid`, and `source /opt/xilinx/xrt/setup.sh`.

2. **Bitstream** (only if `2xor_11bit_hw.xclbin` is lost): run
   `kernel/build_xclbin.sh` on the server (~2 h). It reproduces the exact xclbin
   used from exp40 through exp73.

## Expected output (must match exactly)

```
Unique (acc)              : 109533748     <- must equal KMC exactly
dropped inserts (bin full): 0             <- must be 0
full bins / total bins    : 0 / 32768
occupied slots / capacity : 109533748 / 134217728  (load 81.609%)
Median e2e over 3 runs    : ~10.22 s
```

## Reporting caveats (from SHARED_CONTEXT)

- **CPU baseline must be `-t16`** (the machine is 8P/16T; the host uses all 16
  logical cores + 3 U50s). Older `-t4` comparisons are unfair.
- **Cold bitstream programming** is a one-time ~28.7 s across the 3 cards
  (per bitstream). All benchmarks here were warm; state this when reporting.
- **KMC writes tmp (~4.5 GB) + an output DB to disk; exp73 writes nothing**
  (RAM hash table only). exp73 wins even under this disadvantage.
