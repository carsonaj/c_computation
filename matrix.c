#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "array.h"
#include "polynomial.c"
#include "matrix.h"

#define TRUE 1
#define FALSE 0 


// zero indexed matrix library


//Matrix
// data structure
Matrix mat_create(int rows, int cols) {
    assert(rows*cols <= MAX_SIZE);
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;

    return mat;
}


Matrix mat_zero(int rows, int cols) {
    Matrix mat = mat_create(rows, cols);
    int i;
    for (i=0; i<rows*cols; i++) {
        mat.data[i] = 0.0;
    }

    return mat;
}

double mat_get_element(Matrix mat, int row, int col) {
    return mat.data[row*(mat.cols) + col];
}

void mat_set_element(Matrix *mat, int row, int col, double element) {
    mat->data[row*(mat->cols) + col] = element;
}

Matrix mat_get_rows(Matrix mat, int rows, int *rows_arr) {
    int cols = mat.cols;
    Matrix row_mat = mat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, rows_arr[i], j);
            mat_set_element(&row_mat, i, j, element);
        }
    }

    return row_mat;
}

Matrix mat_get_cols(Matrix mat, int cols, int *cols_arr) {
    int rows = mat.rows;
    Matrix col_mat = mat_create(rows, cols);

    int i, j;
    for (j=0; j<cols; j++) {
        for (i=0; i<rows; i++) {
            double element = mat_get_element(mat, i, cols_arr[j]);
            mat_set_element(&col_mat, i, j, element);
        }
    }

    return col_mat;
}

Matrix mat_join(Matrix A, Matrix B, int axis) {
    // axis = 0 means vertical join
    // axis = 1 means horizontal join
    int m = A.rows;
    int n = A.cols;
    int r = B.rows;
    int s = B.cols;

    Matrix join;
    if (axis==0) {
        assert(n==s);
        join = mat_create(m+r, n);
        int i, j;
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                double elementA = mat_get_element(A, i, j);
                mat_set_element(&join, i, j, elementA);
            }

            for (i=0; i<r; i++) {
                double elementB = mat_get_element(B, i, j);
                mat_set_element(&join, i+m, j, elementB);
            }
        }
    }

    else if (axis==1) {
        assert(m==r);
        join = mat_create(m, n+s);
        int i, j;
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                double elementA = mat_get_element(A, i, j);
                mat_set_element(&join, i, j, elementA);
            }

            for (j=0; j<s; j++) {
                double elementB = mat_get_element(B, i, j);
                mat_set_element(&join, i, j+n, elementB);
            }
        }
    }

    return join;
}

void mat_print(Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, i, j);
            printf("%6.2f ", element);
        }
        printf("\n");
    }
}

// create a copy of a matrix in different memory
Matrix mat_copy(Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    Matrix cpy = mat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, i, j);
            mat_set_element(&cpy, i, j, element);
        }
    }

    return cpy;
}
//------------------------------------------------------------------------

// mathematics

// elementary row operations:

// (type 1) swaps rows i,j
void mat_row_op1(Matrix *mat, int i, int j) {
    int cols = mat->cols;
    int rows_arr[1] = {i};
    Matrix temp_i = mat_get_rows(*mat, 1, rows_arr);
    int col;
    for (col=0; col<cols; col=col+1) {
        double element_j = mat_get_element(*mat, j, col);
        mat_set_element(mat, i, col, element_j);
        double element_i = mat_get_element(temp_i, 0, col);
        mat_set_element(mat, j, col, element_i);
    }
}

// (type 2) multiplies row i by a constant k
void mat_row_op2(Matrix *mat, int i, double k) {
    int cols = mat->cols;
    int j;
    for (j=0; j<cols; j=j+1) {
        double element = k*mat_get_element(*mat, i, j);
        mat_set_element(mat, i, j, element);
    }
}

// (type 3) multiplies row j by constant k and adds it to row i
void mat_row_op3(Matrix *mat, int i, int j, double k) {
    int cols = mat->cols;
    int col;
    for (col=0; col<cols; col=col+1) {
        double element = mat_get_element(*mat, i, col) + k*mat_get_element(*mat, j, col);
        mat_set_element(mat, i, col, element);
    }
}

int mat_equal(Matrix A, Matrix B) {
    int same_rows = (A.rows == B.rows);
    int same_cols = (A.cols == B.cols);
    if (!same_rows || !same_cols) {
        return FALSE;
    }

    int rows = A.rows;
    int cols = A.cols;
    int i;
    for (i=0; i<rows*cols; i=i+1) {
        if (A.data[i] != B.data[i])
            return FALSE;
    }

    return TRUE;
}

// standard matrix product
Matrix mat_product(Matrix A, Matrix B) {
    assert(A.cols == B.rows);
    int m,n,p;
    m = A.rows;
    n = A.cols;
    p = B.cols;

    Matrix mat = mat_create(m, p);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<p; j=j+1) {

            double kron_prod[n];
            int k;
            for (k=0; k<n; k=k+1) {
                kron_prod[k] = mat_get_element(A, i, k)*mat_get_element(B, k, j);
            }

            double element = arr_sum(kron_prod, n);
            mat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

// Hadamard product
Matrix mat_had_product(Matrix A, Matrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);

    int rows = A.rows;
    int cols = A.cols;

    Matrix prod = mat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(A, i, j)*mat_get_element(B, i, j);
            mat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

Matrix mat_scale(double c, Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    Matrix prod = mat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = c*mat_get_element(mat, i, j);
            mat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

Matrix mat_sum(Matrix A, Matrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);
    int m = A.rows;
    int n = A.cols;

    Matrix mat = mat_create(m, n);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<n; j=j+1) {
            double element = mat_get_element(A, i, j) + mat_get_element(B, i, j);
            mat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

Matrix mat_transpose(Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;

    Matrix trans = mat_create(cols, rows);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, i, j);
            mat_set_element(&trans, j, i, element);
        }
    }

    return trans;
}

//----------------------------------------------------------------------

// algorithms

// row echelon form
static void sub_ref(Matrix *mat, int start_row, int start_col) {

    // check if mat is zero matrix
    int rows = mat->rows;
    int cols = mat->cols;
    int val = TRUE;
    double element;
    int row = start_row;
    int col = start_col;
    while (col<cols) {
        row = start_row;
        while (row<rows) {
            element = mat_get_element(*mat, row, col);
            if (element != 0) {
                val = FALSE;
                break;
            }
            row = row+1;
        }
        if (!val)
            break;
        col = col+1;
    }
    if (val)
        return;

    // (row, col) is now pivot;
    // row operations to create a pivot of 1 and zeros below pivot
    mat_row_op2(mat, row, 1.0/element);
    mat_row_op1(mat, start_row, row);

    row = row+1;
    while (row<rows) {
        double k = -mat_get_element(*mat, row, col);
        mat_row_op3(mat, row, start_row, k);
        row = row+1;
    }
    // recursive call on submatrix down-right from pivot
    sub_ref(mat, start_row+1, col+1);
}

void mat_ref(Matrix *mat) {
    sub_ref(mat, 0, 0);
}

// reduced row echelon form
static void sub_rref(Matrix *mat, int start_row, int start_col) {
    // assume matrix is in row echelon form
    // check if mat is zero matrix
    int rows = mat->rows;
    int cols = mat->cols;
    int val = TRUE;
    double element;
    int row = start_row;
    int col = start_col;
    while (col<cols) {
        row = start_row;
        while (row<rows) {
            element = mat_get_element(*mat, row, col);
            if (element != 0) {
                val = FALSE;
                break;
            }
            row = row+1;
        }
        if (!val)
            break;
        col = col+1;
    }
    if (val)
        return;

    // now (row, col) is the pivot
    int prev_row = row-1;
    while (0<=prev_row) {
        element = -mat_get_element(*mat, prev_row, col);
        if (element != 0)
            mat_row_op3(mat, prev_row, row, element);
        prev_row = prev_row-1;
    }

    row = row+1;
    col = col+1;
    sub_rref(mat, row, col);
}

void mat_rref(Matrix *mat) {
    mat_ref(mat);
    sub_rref(mat, 0, 0);
}

// solve system Ax=b
Matrix mat_solve_system(Matrix A, Matrix b) {
    int m = A.rows;
    int n = A.cols;
    int l = b.rows;
    assert(m==n);
    assert(m==l);

    Matrix xsol = mat_join(A, b, 1);
    mat_rref(&xsol);
    int cols_arr[1] = {m};
    Matrix x = mat_get_cols(xsol, 1, cols_arr);

    return x;
}

//-------------------------------------------------------------------









//PolyMatrix
PolyMatrix pymat_create(int rows, int cols) {
    assert(rows*cols <= MAX_SIZE);
    PolyMatrix mat;
    mat.rows = rows;
    mat.cols = cols;

    return mat;
}

PolyMatrix pymat_zero(int rows, int cols) {
    PolyMatrix mat = pymat_create(rows, cols);
    int i;
    for (i=0; i<rows*cols; i++) {
        mat.data[i] = ply_zero();
    }

    return mat;
}

Polynomial pymat_get_element(PolyMatrix mat, int row, int col) {
    return mat.data[row*(mat.cols) + col];
}

void pymat_set_element(PolyMatrix *mat, int row, int col, Polynomial element) {
    mat->data[row*(mat->cols) + col] = element;
}

PolyMatrix pymat_get_rows(PolyMatrix mat, int rows, int *rows_arr) {
    int cols = mat.cols;
    PolyMatrix row_mat = pymat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = pymat_get_element(mat, rows_arr[i], j);
            pymat_set_element(&row_mat, i, j, element);
        }
    }

    return row_mat;
}

PolyMatrix pymat_get_cols(PolyMatrix mat, int cols, int *cols_arr) {
    int rows = mat.rows;
    PolyMatrix col_mat = pymat_create(rows, cols);

    int i, j;
    for (j=0; j<cols; j++) {
        for (i=0; i<rows; i++) {
            Polynomial element = pymat_get_element(mat, i, cols_arr[j]);
            pymat_set_element(&col_mat, i, j, element);
        }
    }

    return col_mat;
}

PolyMatrix pymat_join(PolyMatrix A, PolyMatrix B, int axis) {
    // axis = 0 means vertical join
    // axis = 1 means horizontal join
    int m = A.rows;
    int n = A.cols;
    int r = B.rows;
    int s = B.cols;

    PolyMatrix join;
    if (axis==0) {
        assert(n==s);
        join = pymat_create(m+r, n);
        int i, j;
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Polynomial elementA = pymat_get_element(A, i, j);
                pymat_set_element(&join, i, j, elementA);
            }

            for (i=0; i<r; i++) {
                Polynomial elementB = pymat_get_element(B, i, j);
                pymat_set_element(&join, i+m, j, elementB);
            }
        }
    }

    else if (axis==1) {
        assert(m==r);
        join = pymat_create(m, n+s);
        int i, j;
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                Polynomial elementA = pymat_get_element(A, i, j);
                pymat_set_element(&join, i, j, elementA);
            }

            for (j=0; j<s; j++) {
                Polynomial elementB = pymat_get_element(B, i, j);
                pymat_set_element(&join, i, j+n, elementB);
            }
        }
    }

    return join;
}

PolyMatrix pymat_copy(PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    PolyMatrix cpy = pymat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = pymat_get_element(mat, i, j);
            pymat_set_element(&cpy, i, j, element);
        }
    }

    return cpy;
}

int pymat_equal(PolyMatrix A, PolyMatrix B) {
    int same_rows = (A.rows == B.rows);
    int same_cols = (A.cols == B.cols);
    if (!same_rows || !same_cols) {
        return FALSE;
    }

    int rows = A.rows;
    int cols = A.cols;
    int i;
    for (i=0; i<rows*cols; i=i+1) {
        if (ply_equal(A.data[i], B.data[i])) {
            return FALSE;
        }
    }

    return TRUE;
}

// standard matrix product
PolyMatrix pymat_product(PolyMatrix A, PolyMatrix B) {
    assert(A.cols == B.rows);
    int m,n,p;
    m = A.rows;
    n = A.cols;
    p = B.cols;

    PolyMatrix mat = pymat_create(m, p);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<p; j=j+1) {

            Polynomial kron_prod[n];
            int k;
            for (k=0; k<n; k=k+1) {
                Polynomial Aik = pymat_get_element(A, i, k);
                Polynomial Bkj = pymat_get_element(B, k, j);
                kron_prod[k] = ply_product(Aik, Bkj);
            }

            Polynomial sum[n-1];
            sum[0] = ply_copy(kron_prod[0]);
            int l;
            for (l=1; l<n-1;l++) {
                sum[i] = ply_sum(kron_prod[i-1], kron_prod[i]);
            }

            Polynomial element = sum[n-2];
            pymat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

// Hadamard product
PolyMatrix pymat_had_product(PolyMatrix A, PolyMatrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);

    int rows = A.rows;
    int cols = A.cols;

    PolyMatrix prod = pymat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial Aij = pymat_get_element(A, i, j);
            Polynomial Bij = pymat_get_element(B, i, j);
            Polynomial element = ply_product(Aij, Bij);
            pymat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

PolyMatrix pymat_scale(double c, PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    PolyMatrix prod = pymat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = ply_scale(c, pymat_get_element(mat, i, j));
            pymat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

PolyMatrix pymat_poly_scale(Polynomial p, PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    PolyMatrix prod = pymat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = ply_product(p, pymat_get_element(mat, i, j));
            pymat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

PolyMatrix pymat_sum(PolyMatrix A, PolyMatrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);
    int m = A.rows;
    int n = A.cols;

    PolyMatrix mat = pymat_create(m, n);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<n; j=j+1) {
            Polynomial Aij = pymat_get_element(A, i, j);
            Polynomial Bij = pymat_get_element(B, i, j);
            Polynomial element = ply_sum(Aij, Bij);
            pymat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

PolyMatrix pymat_transpose(PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;

    PolyMatrix trans = pymat_create(cols, rows);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = pymat_get_element(mat, i, j);
            pymat_set_element(&trans, j, i, element);
        }
    }

    return trans;
}
