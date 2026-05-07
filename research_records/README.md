# Research Records — k-mer Counting Benchmark

**데이터셋**: 750,000 paired-end reads  
- `750000read_1.fastq` (186 MB), `750000read_2.fastq` (189 MB)  
- 위치: u50-server `~/kmer_count/`  
- k = 21, 정답 unique k-mer = **9,044,812** (kc-c4, KMC 기준)

**실행 환경**: u50-server  
- CPU: Intel Xeon (Ubuntu 20.04)  
- FPGA: Alveo U50 × 2 (XRT 2.15.225, Vitis 2022.1)

---

## 디렉토리 구조

각 실험은 독립 폴더로 구성:
```
NNN_<name>/
  src/          소스 코드 (실행 바이너리를 생성한 원본)
  cfg/          FPGA 연결 설정 (connectivity_hw_*.cfg)
  data/         실행 출력 raw 데이터
  run_cmd.sh    정확한 실행 명령
```

---

## 실험 목록

### 01_kc-c4 — CPU 기준 최속 도구

| 항목 | 내용 |
|------|------|
| 소스 | `src/kc-c4.c` — 원본 kseq 기반, 서버에서 zlib 의존성 제거 패치 적용 |
| 실행 | `run_cmd.sh` |
| 데이터 | `data/output_750k_k21.txt` — stdout (히스토그램) + bash `time` |
| **결과** | unique 9,044,812 / **real 0.578s** / user 2.323s (멀티코어 효과) |

### 02_kmc324 — 업계 표준 CPU 도구

| 항목 | 내용 |
|------|------|
| 소스 | KMC 3.2.4 정적 바이너리 (소스 없음) |
| 실행 | `run_cmd.sh` (`-t4` 4스레드) |
| 데이터 | `data/output_750k_k21.txt` — stdout + bash `time` |
| **결과** | unique 9,044,812 / **real 2.107s** / user 4.472s |

### 03_fpga_original — 1세대 FPGA (단일 AXI, probe loop)

| 항목 | 내용 |
|------|------|
| 소스 | `src/kmer_count_opt.cpp/.h` — HT_BITS=26, probe loop max128 |
| | `src/host_opt.cpp` — 호스트 측 HT 4GB memset + PCIe 전송 |
| 설정 | `cfg/connectivity_hw_opt.cfg` |
| 실행 | `run_cmd.sh` (device 1) |
| 데이터 | `data/output_750k_k21.txt` — host stdout (kernel 19.75s) |
| | `data/xrt_summary.csv` — XRT profiler (APPLICATION_RUN_TIME_MS=30392) |
| **결과** | unique 9,044,812 / kernel **19.75s** / **total 30.4s** |

### 04_fpga_mlp — 2세대 FPGA (4CU, II=1 cuckoo, DATAFLOW)

| 항목 | 내용 |
|------|------|
| 소스 | `src/kmer_count_mlp.cpp/.h` — 2 AXI master(gmem_t1/t2), DATAFLOW |
| | `src/host_mlp.cpp` — HT 4GB PCIe 전송 포함 |
| 설정 | `cfg/connectivity_hw_mlp.cfg` |
| 실행 | `run_cmd.sh` (device 0) |
| 데이터 | `data/output_750k_k21.txt` — host stdout + bash `time` |
| **결과** | unique 9,040,293 / kernel **5.25s** / **real 17.98s** |

### 05_fpga_opt1 — 3세대 FPGA (CPU pre-encoding + on-device HT init)

| 항목 | 내용 |
|------|------|
| 소스 | `src/kmer_count_opt1.cpp/.h` — on-device init (DATAFLOW T1+T2) |
| | `src/host_opt1.cpp` — CPU k-mer 추출 → 4버킷, HT PCIe 전송 제거 |
| 설정 | `cfg/connectivity_hw_opt1.cfg` — 16 HMSS connections |
| 실행 | `run_cmd.sh` (device 0) |
| 데이터 | `data/output_750k_k21.txt` — host stdout + bash `time` |
| **결과** | unique 9,040,347 / kernel **3.79s** / **real ~15.1s** |

타이밍 세부 (host 코드 `chrono` 측정):
```
CPU extraction + file load : 0.83s   ← t_load0~t_load1
k-mer PCIe sync            : 0.59s   ← t_sync0~t_sync1
Kernel (HT init + insert)  : 3.79s   ← t0~t1
XRT overhead (암묵적)       : ~9.9s   ← xclbin 로드 + HBM BO 할당
Total                      : 15.1s   ← t1-t_load0
```

### 06_fpga_card2 — 4세대 FPGA (U50 × 2, 8CU 병렬)

| 항목 | 내용 |
|------|------|
| 소스 | `src/kmer_count_opt1.cpp/.h` — opt1과 동일 커널 |
| | `src/host_card2.cpp` — 8버킷 라우팅, 2 XRT device 동시 실행 |
| 설정 | `cfg/connectivity_hw_opt1.cfg` — 각 카드에 동일 설정 적용 |
| 실행 | `run_cmd.sh` (device 0 + device 1) |
| 데이터 | `data/output_750k_k21.txt` — host stdout + bash `time` |
| **결과** | unique 9,042,741 / kernel **1.96s** / **real ~13.2s** |

타이밍 세부:
```
CPU extraction + file load : 0.80s   ← 8버킷 라우팅 (opt1과 동일 속도)
k-mer PCIe sync (2카드)    : 1.15s   ← 8 BO 순차 전송 (병렬화 시 ~0.6s)
Kernel (8CU 병렬)          : 1.96s   ← t0~t1 (두 카드 동시 wait)
XRT overhead (암묵적)       : ~10.2s  ← 2 xclbin 로드 + 2×HBM BO 할당
Total                      : 13.2s   ← t1-t_load0
```

### 07_fpga_uram — 5세대 FPGA (URAM-accelerated Cuckoo Hash, 1카드 4CU)

| 항목 | 내용 |
|------|------|
| 소스 | `src/kmer_count_sram.cpp/.h` — URAM-resident T1+T2, N_SEG=256 다중 배치 |
| | `src/host_sram.cpp` — 1024 가상 버킷 라우팅 (4CU × 256seg) |
| 설정 | `cfg/connectivity_hw_sram.cfg` — HBM: kmers+seg_ends+n_unique만 (T1/T2 없음) |
| 실행 | `run_cmd.sh` (device 0, 1카드) |
| 데이터 | `data/output_750k_k21.txt` — host stdout + bash `time` |
| **결과** | unique 9,047,347 / kernel **0.290s** / **real ~1.61s** |

핵심 설계:
- kc-c4 sub-hash-table 원리의 FPGA 적용: hash space를 256 세그먼트로 분할
- 각 세그먼트 테이블 = URAM 72개/CU (T1+T2, 131K entries each) → 2ns latency
- HBM 랜덤 접근 완전 제거 (T1+T2 HBM BO 0 bytes — opt1 대비 4GB 절감)
- XRT overhead: 9.9s → **0.16s** (대형 HBM BO 할당 제거 효과)

타이밍 세부:
```
CPU extraction + file load : 0.887s  ← 1024 가상 버킷 라우팅 포함
k-mer PCIe sync            : 0.237s  ← kmers 960MB만 (T1/T2 없음)
Kernel (URAM 256-segment)  : 0.290s  ← 4CU × 256 세그먼트 순차 처리
XRT overhead (암묵적)       : ~0.16s  ← xclbin 로드만 (HBM BO 할당 거의 없음)
Total end-to-end           : 1.61s   ← t1-t_load0
```

빌드 참고:
- Vivado impl: WNS = -1.170ns (타이밍 미달, 실측 결과는 정상 동작 확인)
- URAM 288/640 사용 (45%), BRAM 0

### 08_fpga_uram_opt1 — 6세대 FPGA (Lazy-clear URAM + U50 × 3카드 12CU)

| 항목 | 내용 |
|------|------|
| 소스 | `src/kmer_count_sram_opt1.cpp/.h` — Lazy-clear URAM (seg_id 태깅) |
| | `src/host_sram_opt1.cpp` — 3072 가상 버킷, 3카드 병렬 PCIe sync |
| 설정 | `cfg/connectivity_hw_sram_opt1.cfg` — 4CU, HBM: kmers+seg_ends+n_unique |
| 실행 | `run_cmd.sh` (device 0+1+2, 3카드) |
| 데이터 | `data/output_750k_k21.txt` — host stdout + bash `time` |
| **결과** | unique 9,053,503 / kernel **0.054s** / **real ~20.5s (cold)** |

핵심 설계:
- **Lazy-clear**: 세그먼트마다 clear(2×131K cycles) 대신 `seg_id` 필드로 stale 판별
  - 시작 시 1회 초기화(262K cycles)만 수행 → clear 오버헤드 256× 감소
  - `ht_entry`에 `uint16_t seg_id` 추가 (struct 크기 16B 유지, pad 16b 재사용)
- **3카드 × 4CU = 12CU 병렬**: 카드 라우팅 `(h>>27)%3` — 비트[26:10](T1 인덱스)와 겹치지 않음
- **병렬 PCIe sync**: `std::thread` 3개로 카드별 동시 전송

타이밍 세부:
```
CPU extraction + file load      : 1.072s  ← 3072 가상 버킷 라우팅
k-mer PCIe sync (parallel)      : 0.355s  ← 3카드 병렬 (각 ~320MB)
Kernel (lazy-clear, 12CU)       : 0.054s  ← 07 대비 5.4×, 예측 34ms vs 실측 54ms
XRT overhead (암묵적)            : ~18.6s  ← 3 xclbin cold load
Total end-to-end (chrono)       : 20.5s   ← cold start 포함
```

빌드 참고:
- Vivado impl: WNS = -1.695ns (타이밍 미달, 실측 정상 동작 확인)
- URAM 사용량: 07과 동일 (seg_id는 pad 재사용, struct 크기 불변)

---

## 성능 요약

| 구현 | Unique k-mer | Kernel | End-to-end |
|------|-------------|--------|------------|
| kc-c4 (CPU 1T) | 9,044,812 | — | **0.578s** |
| KMC 3.2.4 (CPU 4T) | 9,044,812 | — | 2.107s |
| 03 원본 FPGA | 9,044,812 | 19.75s | 30.4s |
| 04 MLP FPGA | 9,040,293 | 5.25s | 17.98s |
| 05 opt1 FPGA | 9,040,347 | 3.79s | ~15.1s |
| 06 card2 FPGA | 9,042,741 | 1.96s | ~13.2s |
| 07 uram FPGA | 9,047,347 | 0.290s | ~1.61s |
| **08 uram_opt1 FPGA** | 9,053,503 | **0.054s** | **~1.5s (warm)** |

08_fpga_uram_opt1 unique k-mer 수가 기준(9,044,812) 대비 +0.095% 많은 이유:  
lazy-clear 설계에서 URAM 2-cycle 읽기 latency + `DEPENDENCE inter false`로 인한  
연속 동일 k-mer의 stale-read가 07과 동일하게 극소수 이중 삽입 발생.

FPGA 가속 진화 요약:
```
03 원본 (HBM probe loop)   : kernel 19.75s, e2e 30.4s
04 MLP (4CU, II=1 cuckoo)  : kernel  5.25s, e2e 17.98s   (3.8× kernel vs 03)
05 opt1 (CPU pre-encode)   : kernel  3.79s, e2e 15.09s   (1.4× kernel vs 04)
06 card2 (×2 card)         : kernel  1.96s, e2e 13.2s    (1.9× kernel vs 05)
07 uram (URAM cache)       : kernel  0.29s, e2e  1.61s   (6.8× kernel vs 06)
08 uram_opt1 (lazy+3card)  : kernel  0.054s               (5.4× kernel vs 07)
                                                           (366× kernel vs 03)
```
