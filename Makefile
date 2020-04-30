all: test

test: polynomial.o counting.o
	gcc -Werror test.c -o test polynomial.o counting.o

polynomial.o: polynomial.c
	gcc -c polynomial.c

counting.o: counting.c
	gcc -c counting.c

matrix.o: matrix.c
	gcc -c matrix.c
