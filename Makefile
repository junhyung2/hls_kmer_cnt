TARGET   ?= sw_emu
PLATFORM  = xilinx_u50_gen3x16_xdma_5_202210_1

KERNEL    = kmer_count
HOST_EXE  = host_$(TARGET)
KRNL_SRC  = src/kmer_count.cpp
HOST_SRC  = src/host.cpp
XO_FILE   = $(KERNEL)_$(TARGET).xo
XCLBIN    = $(KERNEL)_$(TARGET).xclbin

# ── Per-target settings ──────────────────────────────────────────────────────
ifeq ($(TARGET), sw_emu)
  CFG_FILE     = connectivity_sw_emu.cfg
  HT_BITS_HOST = 22
  KRNL_DEFINES = -DSW_EMU
else ifeq ($(TARGET), hw_emu)
  CFG_FILE     = connectivity_hw_emu.cfg
  HT_BITS_HOST = 22
  KRNL_DEFINES = -DHW_EMU
else ifeq ($(TARGET), hw)
  CFG_FILE     = connectivity_hw.cfg
  HT_BITS_HOST = 28
  KRNL_DEFINES =
else
  $(error Unknown TARGET=$(TARGET). Use: sw_emu | hw_emu | hw)
endif

# ── Paths ────────────────────────────────────────────────────────────────────
XILINX_XRT    ?= /opt/xilinx/xrt
XILINX_VITIS  ?= /tools/Xilinx/Vitis/2022.1
VPP            = $(XILINX_VITIS)/bin/v++

VPP_COMMON = --platform $(PLATFORM) --target $(TARGET) --save-temps \
             -I src

# Kernel-side compiler defines passed through v++ HLS front-end
VPP_COMPILE_FLAGS = $(VPP_COMMON) $(KRNL_DEFINES)
VPP_LINK_FLAGS    = $(VPP_COMMON) --config $(CFG_FILE)

# Host compiler flags
CXXFLAGS = -std=c++17 -O2 -Wall \
            -DHT_BITS_HOST=$(HT_BITS_HOST) \
            -I$(XILINX_XRT)/include \
            -I src
LDFLAGS  = -L$(XILINX_XRT)/lib -lxrt_coreutil -lpthread

# ── Rules ────────────────────────────────────────────────────────────────────
.PHONY: all kernel host run clean

all: kernel host

# Step 1: compile kernel source → .xo
kernel: $(XCLBIN)

$(XO_FILE): $(KRNL_SRC) src/kmer_count.h
	@echo ">>> [v++ compile] TARGET=$(TARGET)"
	$(VPP) -c $(VPP_COMPILE_FLAGS) -k $(KERNEL) -o $@ $<

# Step 2: link .xo → .xclbin
$(XCLBIN): $(XO_FILE)
	@echo ">>> [v++ link] TARGET=$(TARGET)"
	$(VPP) -l $(VPP_LINK_FLAGS) -o $@ $<

# Host binary
host: $(HOST_EXE)

$(HOST_EXE): $(HOST_SRC) src/kmer_count.h
	@echo ">>> [g++ host] TARGET=$(TARGET)"
	g++ $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# Emulation config (needed by XRT to find the emulated device)
emconfig.json:
	$(XILINX_VITIS)/bin/emconfigutil --platform $(PLATFORM) --nd 1

# Run (uses small test FASTQs in ../data/)
run: $(HOST_EXE) $(XCLBIN) emconfig.json
	@echo ">>> [run] XCL_EMULATION_MODE=$(TARGET)"
	XCL_EMULATION_MODE=$(TARGET) ./$(HOST_EXE) \
	    $(XCLBIN) \
	    ../data/022075_read1.fastq \
	    ../data/022075_read2.fastq \
	    21

clean:
	rm -f *.xo *.xclbin *.json *.log *.pb *.link_summary *.compile_summary
	rm -rf _x .Xil host_sw_emu host_hw_emu host_hw
