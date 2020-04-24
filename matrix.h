// matrix

typedef struct Matrix_ Matrix;

struct Matrix_ {
    int rows;
    int cols;
    double *data;
};

// structure
Matrix *mat_create_matrix(int rows, int cols);
void mat_delete_matrix(Matrix *mat);
Matrix *mat_zeros(int rows, int cols);
double mat_get_element(Matrix *mat, int row, int col);
void mat_set_element(Matrix *mat, int row, int col, double element);
Matrix *mat_get_rows(Matrix *mat, int rows, int *rows_arr);
Matrix *mat_get_cols(Matrix *mat, int cols, int *cols_arr);
Matrix *mat_join(Matrix *A, Matrix *B, int axis);
void mat_print_matrix(Matrix *mat);
Matrix *mat_copy_matrix(Matrix *mat);

// mathematics
void mat_row_op1(Matrix *mat, int i, int j);
void mat_row_op2(Matrix *mat, int i, double k);
void mat_row_op3(Matrix *mat, int i, int j, double k);
bool mat_equal(Matrix *A, Matrix *B);
Matrix *mat_product(Matrix *A, Matrix *B);
Matrix *mat_had_product(Matrix *A, Matrix *B);
Matrix *mat_scalar_poduct(double c, Matrix *mat);
Matrix *mat_sum(Matrix *A, Matrix *B);
Matrix *mat_transpose(Matrix *mat);

// algorithms
void mat_ref(Matrix *mat);
void mat_rref(Matrix *mat);
Matrix *mat_solve_system(Matrix *A, Matrix *b);
