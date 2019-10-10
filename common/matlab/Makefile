CC = mex
CFLAG = -largeArrayDims

TARGETS = $(basename $(wildcard *.c))

all: $(TARGETS)

loadmatrix: loadmatrix.c io.h
	$(CC) $(CFLAG) loadmatrix.c

loadmatrix_binary: loadmatrix_binary.c io.h binary.h
	$(CC) $(CFLAG) loadmatrix_binary.c

logsum_mex: logsum_mex.c numeric.h
	$(CC) $(CFLAG) logsum_mex.c

savematrix: savematrix.c
	$(CC) $(CFLAG) savematrix.c

get_binary: get_binary.c binary.h
	$(CC) $(CFLAG) get_binary.c

