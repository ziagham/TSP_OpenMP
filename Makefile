CC = gcc

CFLAGS = -O2 -Wall
CFLAGSParallel = -O2 -fopenmp

LIB = -lm

all: StopWatch.o RaceTrap

RaceTrap: RaceTrap.c
	$(CC) $(CFLAGS) RaceTrap.c StopWatch.o -o RaceTrap $(LIB)

RaceTrapParallel: RaceTrapParallel.c
	$(CC) $(CFLAGSParallel) RaceTrapParallel.c StopWatch.o -o RaceTrapParallel $(LIB)

StopWatch.o: StopWatch.c
	$(CC) -c StopWatch.c

clean:
	rm -f *~ *.o core* RaceTrap

cleandata:
	rm -f data/*
