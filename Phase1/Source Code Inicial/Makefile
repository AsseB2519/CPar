CC = gcc
SRC = src/
CFLAGS = -pg

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

gprof:
	$(CC) $(CFLAGSGPROF) $(SRC)MD.cpp -lm -o MD.exe
	./MD.exe < inputdata.txt
	gprof MD.exe gmon.out > gprof.txt
	cat gprof.txt

run:
	./MD.exe < inputdata.txt

clean:
	rm ./MD.exe