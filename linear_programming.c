#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "array.h"
#include "matrix.h"

// A is a mxn matrix, m<=n (m=n has only one single feasible solution,
// so we just worry about the case m<n);
// A corresponds to a linear program in standard form;
// assume rows of A are independent;
// c is a nx1 column matrix, b is an mx1 column matrix,
// but we implement mostly with c transpose, a 1xn row matrix
// to save time/memory

// separate indices for cols that form a basis of column
// space and those that do not by array swapping; this means
// we need m basic indices at the beginning of the array
// since row space and col space have same dimension
int *start_indices(Matrix *A) {
    int m = A->rows;
    int n = A->cols;
    Matrix *cpy = mat_copy_matrix(A);
    mat_rref(cpy);

    int *indices = malloc(n*sizeof(int));
    int l;
    for (l=0; l<n; l++) {
        indices[l] = l;
    }

    int i = 0;
    int j = 0;
    int k = 0;
    while (i<m) {
        while (j<n) {
            if (mat_get_element(cpy, i, j) == 1) {
                arr_intswap(indices, j, k);
                j = j+1;
                k = k+1;
                break;
            }
            else
                j = j+1;
        }
        i = i+1;
    }

    mat_delete_matrix(cpy);
    return indices;
}

Matrix *basic_mat(Matrix *A, int *indices) {
    int m = A->rows;
    int n = A->cols;
    assert(m<n);
    int cols_arr[m];

    int i;
    for (i=0; i<m; i++) {
        cols_arr[i] = indices[i];
    }

        Matrix *B = mat_get_cols(A, m, cols_arr);
        return B;
}

Matrix *nonbasic_mat(Matrix *A, int *indices) {
    int m = A->rows;
    int n = A->cols;
    assert(m<n);
    int cols_arr[n-m];

    int i;
    for(i=0; i<n-m; i++) {
        cols_arr[i] = indices[m+i];
    }

    Matrix *N = mat_get_cols(A, n-m, cols_arr);
    return N;
}

Matrix *basic_cost_trans(Matrix *c_trans, int m, int *indices) {
    int n = c_trans->cols;
    assert(m<n);
    int cols_arr[m];

    int i;
    for(i=0; i<m; i++) {
        cols_arr[i] = indices[i];
    }

    Matrix *cb_trans = mat_get_cols(c_trans, m, cols_arr);
    return cb_trans;
}

Matrix *nonbasic_cost_trans(Matrix *c_trans, int m, int *indices) {
    int n = c_trans->cols;
    assert(m<n);
    int cols_arr[n-m];

    int i;
    for (i=0; i<n-m; i++) {
        cols_arr[i] = indices[m+i];
    }

    Matrix *cn_trans = mat_get_cols(c_trans, n-m, cols_arr);
    return cn_trans;
}

int entering_col(Matrix *B, Matrix *N, Matrix *cb_trans, Matrix *cn_trans) {
    int m = B->rows;
    int n = (N->cols)+m;

    Matrix *cb = mat_transpose(cb_trans);
    Matrix *B_trans = mat_transpose(B);
    Matrix *y = mat_solve_system(B_trans, cb);
    Matrix *y_trans = mat_transpose(y);

    mat_delete_matrix(cb);
    mat_delete_matrix(B_trans);
    mat_delete_matrix(y);

    Matrix *y_transN = mat_product(y_trans, N);
    Matrix *minus_cn_trans = mat_scalar_poduct(-1.0, cn_trans);
    Matrix *zn_trans = mat_sum(y_transN, minus_cn_trans);

    mat_delete_matrix(y_trans);
    mat_delete_matrix(y_transN);
    mat_delete_matrix(minus_cn_trans);

    int val = 0;
    int i;
    i = 0;
    while (i<n-m) {
        if (mat_get_element(zn_trans, 0, i) < 0) {
            val = 1;
            break;
        }
        i = i+1;
    }

    // the index of the associated column in A is indices[m+i]
    if (val==0)
        return -1;
    else
        return i;
}

// returns leaving column and updates xb
int leaving_col(Matrix *B, Matrix *N, Matrix *xb, int e) {
    int m = B->rows;
    int cols_arr[1] = {e};
    Matrix *Ae = mat_get_cols(N, 1, cols_arr);
    Matrix *d = mat_solve_system(B, Ae);

    double t;
    double vals[m];
    int i;
    int j;
    for (i=0; i<m; i++) {
        double elm_xb = mat_get_element(xb, i, 0);
        double elm_d = mat_get_element(d, i, 0);
        if (elm_xb == 0) {
            t = 0.0;
            mat_set_element(xb, i, 0, t);
            return i;
        }

        if (elm_d != 0)
            vals[i] = elm_xb/elm_d;
        else
            vals[i] = -1.0;
    }

    double max = vals[0];
    for (i=0; i<m; i++) {
        if (vals[i] >= max)
            max = vals[i];
    }
    if (max < 0)
        return -1;

    else {
        double min = max;
        j = 0;
        for (i=0; i<m; i++) {
            if (vals[i]<=min && vals[i]>=0)
                min = vals[i];
                j = i;
        }
    }

    t = vals[j];
    Matrix *minus_td = mat_scalar_poduct(-t, d);
    Matrix *xb_t = mat_sum(xb, minus_td);

    mat_delete_matrix(d);
    mat_delete_matrix(minus_td);

    mat_set_element(xb_t, j, 0, t);
    for (i=0; i<m; i++) {
        double element = mat_get_element(xb_t, i, 0);
        mat_set_element(xb, i, 0, element);
    }

    mat_delete_matrix(xb_t);
    return j;
}

void update(Matrix *B, Matrix *N, Matrix *cb_trans, Matrix *cn_trans, int l, int e) {
    int m = B->rows;
    double temp1[m];
    int i;
    for (i=0; i<m; i++) {
        temp1[i] = mat_get_element(B, i, l);
    }

    for (i=0; i<m; i++) {
        double element = mat_get_element(N, i, e);
        mat_set_element(B, i, l, element);
        mat_set_element(N, i, e, temp1[i]);
    }

    double temp2 = mat_get_element(cb_trans, 0, l);
    double element = mat_get_element(cn_trans, 0, e);
    mat_set_element(cb_trans, 0, l, element);
    mat_set_element(cn_trans, 0, e, temp2);
}

int alg(Matrix *B, Matrix *N, Matrix *cb_trans, Matrix *cn_trans, Matrix *xb, int *indices) {
    int m = B->rows;
    int e = entering_col(B, N, cb_trans, cn_trans);
    if (e==-1)
        return 1;
    assert(e>=0);
    int l = leaving_col(B, N, xb, indices, e);
    if (l== -1)
        return -1;

    arr_intswap(indices, l, m+e);
    update(B, N, cb_trans, cn_trans, l, e);
    alg(B, N, cb_trans, cn_trans, xb, indices);
}

Matrix *lp_simplex_method(Matrix *c, Matrix *A, Matrix *b) {
    int m = A->rows;
    int n = A->cols;
    assert(m<=n);
    if (m==n) {
        Matrix *x = mat_solve_system(A, b);
        return x;
    }

    else {
        int *indices = start_indices(A);
        Matrix *B = basic_mat(A, indices);
        Matrix *N = nonbasic_mat(A, indices);

        Matrix *c_trans = mat_transpose(c);
        Matrix *cb_trans = basic_cost_trans(c_trans, m, indices);
        Matrix *cn_trans = nonbasic_cost_trans(c_trans, m, indices);
        mat_delete_matrix(c_trans);
        Matrix *xb = mat_solve_system(B, b);
        int val = alg(B, N, cb_trans, cn_trans, xb, indices);
        assert(val==1);
        Matrix *x = mat_zeros(n, 1);
        int i;
        for (i=0; i<m; i++) {
            double element = mat_get_element(xb, i, 0);
            mat_set_element(x, indices[i], 0, element);
        }

        mat_delete_matrix(B);
        mat_delete_matrix(N);
        mat_delete_matrix(cb_trans);
        mat_delete_matrix(cn_trans);
        mat_delete_matrix(xb);
        free(indices);
        return x;
    }
}
