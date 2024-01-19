# CC1 = gcc
# CC2 = nvcc
# SRC = src/
# CFLAGS = -g -Wall -pg -O2 -ftree-vectorize -mfpmath=sse -march=x86-64 -fno-omit-frame-pointer
# CUDAFLAGS = -O2 -g -std=c++11 -arch=sm_35 -Wno-deprecated-gpu-targets # -G -allow-unsupported-compiler -arch=sm_35

# .DEFAULT_GOAL = all

# all: MDseq.exe MDCuda.exe

# MDseq.exe: $(SRC)/MDseq.cpp
# 	module load gcc/11.2.0;\
# 	$(CC1) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

# MDCuda.exe: $(SRC)/MDCuda.cu
# 	module load cuda/11.3.1;\
# 	module load gcc/9.3.0;\
# 	$(CC2) $(CUDAFLAGS) $(SRC)MDCuda.cu -lm -o MDCuda.exe

# clean:
# 	find . -type f \( ! -path "./src/*" ! -name "Makefile" !  -name "inputdata.txt" \) -exec rm -v {} \;

# runseq: MDseq.exe
# 	./MDseq.exe < inputdata.txt

# runcuda: MDCuda.exe
# 	./MDCuda.exe < inputdata.txt

################################################################################
# Makefile for general code snippets
#
# by André Pereira
################################################################################

SHELL = /bin/sh
BIN_NAME = cuda

CXX = nvcc
LD  = nvcc

CXXFLAGS   = -O2 -g -std=c++11 -arch=sm_35 -Wno-deprecated-gpu-targets 

SRC_DIR = src
BIN_DIR = bin
BUILD_DIR = build
SRC = $(wildcard $(SRC_DIR)/*.cu)
OBJ = $(patsubst src/%.cu,build/%.o,$(SRC))
BIN = $(BIN_NAME)

vpath %.cu $(SRC_DIR)

################################################################################
# Rules
################################################################################

.DEFAULT_GOAL = all

$(BUILD_DIR)/%.o: %.cu
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

$(BIN_DIR)/$(BIN_NAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJ) $(LIBS)

checkdirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

all: checkdirs $(BIN_DIR)/$(BIN_NAME)

# load:
# 	module load gcc/7.2.0;
# 	module load cuda/11.3.1;

run:
	sbatch --partition cpar --constraint=k20 --ntasks=16 --time=5:00 ./run.sh
# sbatch run.sh

clean:
	find . -type f \( ! -path "./src/*" ! -name "Makefile" ! -name "run.sh" ! -name "inputdata.txt" \) -exec rm -v {} \;
