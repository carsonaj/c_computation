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
    Polynomial f = ply_create(2);
    f.coefs[0] =6.0;
    f.coefs[1] =7.0;
    f.coefs[2] =1.0;


    Polynomial g = ply_create(2);
    g.coefs[0] = -6.0;
    g.coefs[1] = -5.0;
    g.coefs[2] = 1.0;

    Polynomial h = ply_create(1);
    h.coefs[0] = -1.0;
    h.coefs[1] = -1.0;

    Polynomial k = ply_create(1);
    k.coefs[0] = 2.0;
    k.coefs[1] = 2.0;


    PolyMatrix A = pymat_create(1,2);
    PolyMatrix B = pymat_create(1,2);

    pymat_set_element(A, 0, 0, f);
    pymat_set_element(A, 0, 1, g);
    pymat_set_element(B, 0, 0, h);
    pymat_set_element(B, 0, 1, k);

    PolyMatrix sum = pymat_sum(A,B);

    int i;
    for (i=0; i<2; i++) {
        ply_print(pymat_get_element(sum,0,i));
    }

    for (i=0; i< 2; i++) {
        ply_delete(sum.data[i]);
    }

    pymat_delete(sum);

    ply_delete(f);
    ply_delete(g);
    ply_delete(h);
    ply_delete(k);

    //PROBLEM see ply_gcd
    //Polynomial *res = ply_gcd(f,g); //PROBLEM

    //printf("gcd is \n");
    //ply_print(res[0]);

    //printf("res[1].deg is %d", res[1].deg);


    //printf("Bezout coefs are \n");
    //ply_print(res[1]);

    //printf("\n and \n" );
    //ply_print(res[2]);
/*
    ply_delete(f);
    ply_delete(g);

    int i;
    for (i=0; i<2; i++) {
        ply_delete(res[i]);
    }*/


    return 0;
}
