// polynomial

typedef struct Polynomial Polynomial;

struct Polynomial {
    int deg;
    double *coefs;
};

// structure
Polynomial ply_create(int deg);
void ply_delete(Polynomial poly);
int ply_get_deg(Polynomial poly);
double ply_get_coef(Polynomial poly, int i);
void ply_set_coef(Polynomial poly, int i, double val);
Polynomial ply_copy(Polynomial p);
void ply_print(Polynomial poly);

// mathematics

// algebra
Polynomial ply_zero();
int ply_is_zero(Polynomial p);
Polynomial ply_sum(Polynomial poly1, Polynomial poly2);
Polynomial ply_product(Polynomial poly1, Polynomial poly2);
Polynomial ply_scale(double s, Polynomial p);
Polynomial *ply_division(Polynomial f, Polynomial g);
Polynomial *ply_gcd(Polynomial p, Polynomial q);
Polynomial *ply_modulo(Polynomial *p, Polynomial *modulus);

// analysis
double ply_evaluate(Polynomial poly, double x);
Polynomial ply_differentiate(Polynomial poly, int n);

// families of polynomials
Polynomial ply_monomial(int n);
Polynomial ply_legendre(int n);
