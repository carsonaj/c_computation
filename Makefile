all: test

test: polynomial.o counting.o matrix.o
	gcc -Wall -Werror test.c -o test polynomial.o counting.o -g -std=c99

matrix.o: matrix.c polynomial.h
	gcc -c matrix.c

polynomial.o: polynomial.c
	gcc -c polynomial.c

counting.o: counting.c
	gcc -c counting.c
