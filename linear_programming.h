Matrix *lp_simplex_method(Matrix *c, Matrix *A, Matrix *b);
int *start_indices(Matrix *A);
Matrix *basic_mat(Matrix *A, int *indices);
Matrix *nonbasic_mat(Matrix *A, int *indices);
Matrix *basic_cost_trans(Matrix *c_trans, int m, int *indices);
Matrix *nonbasic_cost_trans(Matrix *c_trans, int m, int *indices);
int entering_col(Matrix *B, Matrix *N, Matrix *cb_trans, Matrix *cn_trans);
int leaving_col(Matrix *B, Matrix *N, Matrix *xb, int *indices, int e);
