CC = gcc
SRC = src/
CFLAGS = -g -Wall -pg -O2 -ftree-vectorize -mfpmath=sse -march=x86-64 -fno-omit-frame-pointer

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	find . -type f \( ! -path "./src/*" ! -name "Makefile" !  -name "inputdata.txt" ! -name "test.sh" \) -exec rm -v {} \;

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=32;
	./MDpar.exe < inputdata.txt
