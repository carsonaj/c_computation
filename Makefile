all: test

test: polynomial.o counting.o matrix.o
	gcc -Wall -Werror test.c -o test counting.o array.o matrix.o -g -std=c99

matrix.o: matrix.c
	gcc -c matrix.c polynomial.o

polynomial.o: polynomial.c
	gcc -c polynomial.c

counting.o: counting.c
	gcc -c counting.c

array.o: array.c
	gcc -c array.c
