all: test

test: polynomial.o counting.o
	gcc test.c $^ -o test -lm

polynomial.o: counting.o polynomial.c polynomial.h
	gcc -c polynomial.c $< -o polynomial.o

counting.o: counting.c counting.h
	gcc -c $< -o counting.o

be_poly: be_counting
	gcc -c polynomial.c -o polynomial.o
	ls
	gcc test.c polynomial.o counting.o -lm -o test

be_counting:
	gcc -c counting.c -o counting.o
