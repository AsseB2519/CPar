CC = gcc
SRC = src/
CFLAGS1 = -ftree-vectorize -O2 -msse4 -mavx
CFLAGS2 = -ftree-vectorize -O3 -msse4 -mavx
CFLAGSGPROF = -pg -O2 -ftree-vectorize -msse4

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS1) $(SRC)MD.cpp -lm -o MD.exe

srun1:
	$(CC) $(CFLAGS1) $(SRC)MD.cpp -lm -o MD.exe
	srun --partition=cpar perf stat -e cache-references,cache-misses -M cpi ./MD.exe < inputdata.txt

srun2:
	$(CC) $(CFLAGS2) $(SRC)MD.cpp -lm -o MD.exe
	srun --partition=cpar perf stat -e cache-references,cache-misses -M cpi ./MD.exe < inputdata.txt

gprof:
	$(CC) $(CFLAGSGPROF) $(SRC)MD.cpp -lm -o MD.exe
	./MD.exe < inputdata.txt
	gprof MD.exe gmon.out > gprof.txt
	cat gprof.txt

run:
	./MD.exe < inputdata.txt

clean:
	find . -type f \( ! -path "./src/*" ! -name "Makefile" ! -name "inputdata.txt" \) -exec rm -v {} \;