all: tester

tester: polynomial.o counting.o matrix.o number_field.o
	gcc -Wall -Werror tester.c -o tester matrix.o array.o counting.o number_field.o -g -std=c99

number_field.o: number_field.c  
	gcc -c number_field.c polynomial.o

matrix.o: matrix.c
	gcc -c matrix.c polynomial.o

polynomial.o: polynomial.c
	gcc -c polynomial.c

counting.o: counting.c
	gcc -c counting.c

array.o: array.c
	gcc -c array.c
