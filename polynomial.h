#include "matrix.h"

// polynomial
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#define MAX_DEG 100

typedef struct Polynomial Polynomial;

struct Polynomial {
    int deg;
    double coefs[MAX_DEG+1];
};

// structure
Polynomial ply_create(int deg);
Polynomial ply_create_coef(int deg, double *coefs);
int ply_get_deg(Polynomial poly);
double ply_get_coef(Polynomial poly, int i);
void ply_set_coef(Polynomial *poly, int i, double val);
Polynomial ply_copy(Polynomial p);
void ply_print(Polynomial poly);

// mathematics

// algebra
int ply_equal(Polynomial p, Polynomial q);
Polynomial ply_zero();
int ply_is_zero(Polynomial p);
Polynomial ply_sum(Polynomial poly1, Polynomial poly2);
Polynomial ply_product(Polynomial poly1, Polynomial poly2);
Polynomial ply_scale(double s, Polynomial p);
Polynomial ply_neg(Polynomial p);
PolyMatrix ply_division(Polynomial f, Polynomial g);
PolyMatrix ply_gcd(Polynomial p, Polynomial q);
Polynomial pymod_reduce(Polynomial p, Polynomial m);
Polynomial pymod_inv(Polynomial p, Polynomial m);
Polynomial pymod_sum(Polynomial p, Polynomial q, Polynomial m);
Polynomial pymod_product(Polynomial p, Polynomial q, Polynomial m);

// analysis
double ply_evaluate(Polynomial poly, double x);
Polynomial ply_differentiate(Polynomial poly, int n);

// families of polynomials
Polynomial ply_monomial(int n);
Polynomial ply_legendre(int n);

#endif
