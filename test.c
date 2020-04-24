#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "matrix.h"
#include "linear_programming.h"
#include "polynomial.h"
#include "counting.h"

int main() {

    // LinearProgramming Test
/*
    Matrix *A = mat_create_matrix(3, 6);
    double dats1[18] = {1.0,1.0,2.0,1.0,0.0,0.0,  2.0,0.0,3.0,0.0,1.0,0.0,  2.0,1.0,3.0,0.0,0.0,1.0};
    int i;
    for (i=0; i<18; i++) {
        A->data[i] = dats1[i];
    }
    int m = A->rows;
    int n = A->cols;

    Matrix *b = mat_create_matrix(3, 1);
    double dats2[3] = {4.0,5.0,7.0};
    for (i=0; i<3; i++) {
        b->data[i] = dats2[i];
    }

    Matrix *c = mat_create_matrix(6, 1);
    double dats3[6] = {3.0,2.0,4.0,0.0,0.0,0.0};
    for (i=0; i<6; i++) {
        c->data[i] = dats3[i];
    }


    Matrix *x = lp_simplex_method(c, A, b);
    mat_print_matrix(x);

    //int indices[6] = {3,4,5,0,1,2};


    Matrix *c_trans = mat_transpose(c);
    mat_delete_matrix(c);



  int *indices = start_indices(A);
    for (i=0; i<6; i++) {
        printf("%d ", indices[i]);
    }
    printf("\n");



    Matrix *B = basic_mat(A, indices);
    mat_print_matrix(B);
    printf("\n");

    Matrix *N = nonbasic_mat(A, indices);
    mat_print_matrix(N);
    printf("\n");

    Matrix *cb_trans = basic_cost_trans(c_trans, m, indices);
    Matrix *cn_trans = nonbasic_cost_trans(c_trans, m, indices);

    mat_print_matrix(cb_trans);
    printf("\n");
    mat_print_matrix(cn_trans);
    printf("\n");

    Matrix *xb = mat_solve_system(B, b);
    mat_print_matrix(xb);
    printf("\n");




    int e = entering_col(B, N, cb_trans, cn_trans);
    printf("e is %d\n", e);
    printf("\n");


    int l = leaving_col(B, N, xb, indices, e);
    printf("%d\n", l);
    printf("\n");




    m = B->rows;
    int cols_arr[1] = {e};
    Matrix *Ae = mat_get_cols(N, 1, cols_arr);
    mat_print_matrix(Ae);
    printf("\n" );
    Matrix *d = mat_solve_system(B, Ae);
    mat_print_matrix(d);
*/
//----------------------------------------------------------------------


/*    // Matrix Test

    Matrix *B = mat_create_matrix(3, 3, false);
    double dats2[9] = {1,2,3,  2,-1,1,  3,0,-1};
    int i;
    for (i=0; i<9; i++) {
        B->data[i] = dats2[i];
    }

    printf("B is \n");
    mat_print_matrix(B);


    Matrix *b = mat_create_matrix(3, 1, false);
    double dats3[3] = {9,8,3};
    for (i=0; i<3; i++) {
        b->data[i] = dats3[i];
    }



    Matrix *x = mat_solve_system(B, b);
    printf("\n" );
    mat_print_matrix(x);
*/

//-----------------------------------------------------------------


    // polynomial test

/*
    double coef1[5] = {0,1,2, 2,5};
    double coef2[9] = {0,0,0, 10,5,4, 0,0,1};

    Polynomial *p1 = ply_create_poly(4);

    int i;
    for (i=0;i<=4;i++) {
        ply_set_coef(p1, i, coef1[i]);
    }

    Polynomial *p2 = ply_scale(1.0/2, p1);

    ply_print_poly(p1);
    ply_print_poly(p2);
*/

    //(1.0/pow(2.0, n))*(1.0/cnt_factorial(n, 1))*pow(-1.0, exp)*cnt_combination(n, k);


    Polynomial *leg50 = ply_legendre(50);
    printf("\n\n");

    ply_print_poly(leg50);


    return 0;

}
