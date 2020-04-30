#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "matrix.h"
#include "polynomial.h"


int main() {
/*
    Polynomial *z = ply_zero();
    Polynomial *o = ply_monomial(0);

    Polynomial *sum = ply_sum(o,z);
    ply_print_poly(sum);

*/
    Polynomial *f = ply_create_poly(2);
    f->coefs[0] =6;
    f->coefs[1] =7;
    f->coefs[2] =1;


    Polynomial *g = ply_create_poly(2);
    g->coefs[0] = -6;
    g->coefs[1] = -5;
    g->coefs[2] = 1;


    Polynomial **res = ply_gcd(f,g);

    printf("gcd is \n");
    ply_print_poly(res[0]);

    printf("Bezout coefs are \n");
    ply_print_poly(res[1]);

    printf("\n and \n" );
    ply_print_poly(res[2]);


    return 0;
}
