// polynomial

typedef struct Polynomial Polynomial;

struct Polynomial {
    int deg;
    double *coefs;
};

// structure
Polynomial *ply_create_poly(int deg);
void ply_delete_poly(Polynomial *poly);
double ply_get_coef(Polynomial *poly, int i);
void ply_set_coef(Polynomial *poly, int i, double val);
void ply_print_poly(Polynomial *poly);

// mathematics
double ply_evaluate(Polynomial *poly, double x);
Polynomial *ply_zero();
Polynomial *ply_sum(Polynomial *poly1, Polynomial *poly2);
Polynomial *ply_product(Polynomial *poly1, Polynomial *poly2);
Polynomial *ply_scale(double s, Polynomial *p);
Polynomial *ply_differentiate(Polynomial *poly, int n);

// families of polynomials
Polynomial *ply_monomial(int n);
Polynomial *ply_legendre(int n);
