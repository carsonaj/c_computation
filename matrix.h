#ifndef MATRIX_H
#define MATRIX_H
#define MAX_SIZE 500

typedef struct Matrix Matrix;
typedef struct PolyMatrix PolyMatrix;

#include "polynomial.h"

struct Matrix {
    int rows;
    int cols;
    double data[MAX_SIZE];
};

struct PolyMatrix {
    int rows;
    int cols;
    Polynomial data[MAX_SIZE];
};

// Matrix
// structure
Matrix mat_create(int rows, int cols);
Matrix mat_zero(int rows, int cols);
double mat_get_element(Matrix mat, int row, int col);
void mat_set_element(Matrix *mat, int row, int col, double element);
Matrix mat_get_rows(Matrix mat, int rows, int *rows_arr);
Matrix mat_get_cols(Matrix mat, int cols, int *cols_arr);
Matrix mat_join(Matrix A, Matrix B, int axis);
void mat_print(Matrix mat);
Matrix mat_copy(Matrix mat);

// mathematics
void mat_row_op1(Matrix *mat, int i, int j);
void mat_row_op2(Matrix *mat, int i, double k);
void mat_row_op3(Matrix *mat, int i, int j, double k);
int mat_equal(Matrix A, Matrix B);
Matrix mat_product(Matrix A, Matrix B);
Matrix mat_had_product(Matrix A, Matrix B);
Matrix mat_scale(double c, Matrix mat);
Matrix mat_sum(Matrix A, Matrix B);
Matrix mat_transpose(Matrix mat);

// algorithms
void mat_ref(Matrix *mat);
void mat_rref(Matrix *mat);
Matrix mat_solve_system(Matrix A, Matrix b);


//------------------------------------------------------------

// PolyMatrix
//structure
PolyMatrix pymat_create(int rows, int cols);
PolyMatrix pymat_zero(int rows, int cols);
Polynomial pymat_get_element(PolyMatrix mat, int row, int col);
void pymat_set_element(PolyMatrix *mat, int row, int col, Polynomial element);
PolyMatrix pymat_get_rows(PolyMatrix mat, int rows, int *rows_arr);
PolyMatrix pymat_get_cols(PolyMatrix mat, int cols, int *cols_arr);
PolyMatrix pymat_join(PolyMatrix A, PolyMatrix B, int axis);
//void pymat_print(PolyMatrix mat);
PolyMatrix pymat_copy(PolyMatrix mat);

// mathematics
int pymat_equal(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_product(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_had_product(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_scale(double c, PolyMatrix mat);
PolyMatrix pymat_poly_scale(Polynomial p, PolyMatrix mat);
PolyMatrix pymat_sum(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_transpose(PolyMatrix mat);

#endif
