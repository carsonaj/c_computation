#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "polynomial.h"
#include "number_field.h"


int main() {


    Polynomial m = ply_create(2);
    m.coefs[0] =-2.0;
    m.coefs[1] =0.0;
    m.coefs[2] =1.0;

    //Polynomial f = ply_zero();

    Polynomial p = ply_create(1);
    p.coefs[0] = 1.0;
    p.coefs[1] = 2.0;

    NFNumber x;
    x.min_poly = m;
    x.number = p;

    printf("x is \n" );
    nf_print(x);
    printf("x inverse is \n" );
    nf_print(nf_inv(x));

    //PolyMatrix pair = ply_division(f,g);
    //ply_print(pymat_get_element(pair,0,0));
    //ply_print(pymat_get_element(pair,0,1));

    /*Polynomial h = ply_create(3);
    h.coefs[0] = 1.0;
    h.coefs[1] = 1.0;
    h.coefs[2] = 1.0;
    h.coefs[3] = 1.0;

    Polynomial k = ply_create(3);
    k.coefs[0] = 1.0;
    k.coefs[1] = 0.0;
    k.coefs[2] = 0.0;
    k.coefs[3] = 1.0;*/

    /*printf("f is \n" );
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






    /*//PROBLEM see ply_gcd
    PolyMatrix res = ply_gcd(p,m); //PROBLEM
    printf("gcd is \n");
    ply_print(pymat_get_element(res, 0, 0));


    Polynomial res01 = pymat_get_element(res,0,1);
    Polynomial res02 = pymat_get_element(res,0,2);


    printf("Bezout coefs are \n");
    ply_print(res01);

    printf("\n and \n" );
    ply_print(res02);

    printf("pmodm is  is \n" );
    ply_print(pymod_reduce(p,m));

    printf("inverse is \n" );
    ply_print(pymod_inv(p,m));*/



    return 0;
}
