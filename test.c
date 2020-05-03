#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "matrix.h"
#include "polynomial.h"


int main() {


    Polynomial f = ply_create(3);
    f.coefs[0] =-8.0;
    f.coefs[1] =0.0;
    f.coefs[2] =0.0;
    f.coefs[3] =1.0;

    //Polynomial f = ply_zero();

    Polynomial g = ply_create(1);
    g.coefs[0] = -2.0;
    g.coefs[1] = 1.0;

    //PolyMatrix pair = ply_division(f,g);
    //ply_print(pymat_get_element(pair,0,0));
    //ply_print(pymat_get_element(pair,0,1));
/*
    Polynomial h = ply_create(1);
    h.coefs[0] = 1.0;
    h.coefs[1] = 1.0;

    Polynomial k = ply_create(1);
    k.coefs[0] = 2.0;
    k.coefs[1] = 2.0;

    printf("f is \n" );
    ply_print(f);
    printf("g is \n" );
    ply_print(g);
    printf("h is \n" );
    ply_print(h);
    printf("k is \n" );
    ply_print(k);

    printf("f+g is \n" );
    Polynomial sum = ply_sum(f,g);
    ply_print(sum);

    printf("the sum is\n" );

    ply_print(ply_sum(h,k));
    printf("the prod is\n" );
    ply_print(ply_product(h,k));


    PolyMatrix A = pymat_create(1,2);
    PolyMatrix B = pymat_create(1,2);

    pymat_set_element(A, 0, 0, f);
    pymat_set_element(A, 0, 1, g);
    pymat_set_element(B, 0, 0, h);
    pymat_set_element(B, 0, 1, k);

//    PolyMatrix sum = pymat_sum(A,B);

    int i;
    for (i=0; i<2; i++) {
        ply_print(pymat_get_element(sum,0,i);
    }

    */






    //PROBLEM see ply_gcd
    printf("line 68\n" );
    PolyMatrix res = ply_gcd(f,g); //PROBLEM
    printf("line 69\n" );
    printf("gcd is \n");
    ply_print(pymat_get_element(res, 0, 0));


    Polynomial res01 = pymat_get_element(res,0,1);
    Polynomial res02 = pymat_get_element(res,0,2);
    printf("res[0][1].deg is %d\n", res01.deg);



    printf("Bezout coefs are \n");
    ply_print(res01);

    printf("\n and \n" );
    ply_print(res02);

    

    return 0;
}
