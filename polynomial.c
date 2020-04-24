#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "counting.h"
#include "polynomial.h"

// helpful functions
static Polynomial *max_deg(Polynomial *p1, Polynomial *p2) {
    int deg1 = p1->deg;
    int deg2 = p2->deg;
    if (deg1>deg2)
        return p1;
    else
        return p2;
}

static Polynomial *min_deg(Polynomial *p1, Polynomial *p2) {
    int deg1 = p1->deg;
    int deg2 = p2->deg;
    if (deg1>deg2)
        return p2;
    else
        return p1;
}

// data structure
Polynomial *ply_create_poly(int deg) {
    Polynomial *poly = malloc(sizeof(Polynomial));
    poly->deg = deg;

    double *coef = malloc((deg+1)*sizeof(double));
    poly->coefs = coef;

    return poly;
}

void ply_delete_poly(Polynomial *poly) {
    free(poly->coefs);
    free(poly);
}

double ply_get_coef(Polynomial *poly, int i) {
    return poly->coefs[i];
}

void ply_set_coef(Polynomial *poly, int i, double val) {
    poly->coefs[i] = val;
}

void ply_print_poly(Polynomial *poly) {
    int deg = poly->deg;
    double coef0 = ply_get_coef(poly, 0);
    double coefn = ply_get_coef(poly, deg);
    if (deg!=0) {
        assert(coefn!=0);
        int k;
        int i;
        for (i=0; i<=deg; i++) {
            double coef = ply_get_coef(poly, i);
            if (coef != 0) {
                k = i;
                break;
            }
        }

        if (k<deg) {
            double coef = ply_get_coef(poly, k);
            printf("%.2fx^%d ", coef, k);
            for (i=k+1; i<=deg-1; i++) {
                double coef = ply_get_coef(poly, i);
                if (coef>0) {
                    printf("+ %.2fx^%d ", coef, i);
                }
                else if (coef<0)
                    printf("- %.2fx^%d ", -coef, i);
                }
            if (coefn>0)
                printf("+ %.2fx^%d\n", coefn, deg);
            else
                printf("- %.2fx^%d\n", -coefn, deg);
        }

        else if (k==deg) {
            printf("%.3fx^%d\n", coefn, deg);
        }
    }

    else if (deg==0) {
        printf("%.2f\n", coef0);
    }
}
//----------------------------------------------------------------------------

//mathematics

// Horner's method
double ply_evaluate(Polynomial *poly, double x) {
    int deg = poly->deg;
    double coefn = ply_get_coef(poly, deg);
    assert(coefn!=0);

    double val = coefn;
    if (deg==0) {;
        return val;
    }
    else if (deg>0) {
        while (deg>0) {
            val = val*x + ply_get_coef(poly, deg-1);
            deg = deg-1;
        }
    }

    return val;
}

// algebra

Polynomial *ply_zero() {
    Polynomial *z = ply_create_poly(0);
    z->coefs[0] = 0.0;
    return z;
    }

Polynomial *ply_sum(Polynomial *poly1, Polynomial *poly2) {
    Polynomial *p1 = max_deg(poly1, poly2);
    Polynomial *p2 = min_deg(poly1, poly2);
    int deg1 = p1->deg;
    int deg2 = p2->deg;

    double sum_coefs[deg1+1];
    int i;
    for (i=0; i<=deg2; i++) {
        sum_coefs[i] = ply_get_coef(poly1, i) + ply_get_coef(poly2, i);
    }
    for (i=deg2+1; i<=deg1; i++) {
        sum_coefs[i] = ply_get_coef(p1, i);
    }
    int deg = deg1;
    while (deg>=0) {
        if (sum_coefs[deg]!=0)
            break;
        if (deg==0)
            break;
        deg = deg-1;
    }

    Polynomial *sum_poly = ply_create_poly(deg);
    for (i=0; i<=deg; i++) {
        double coef = sum_coefs[i];
        ply_set_coef(sum_poly, i, coef);
    }

    return sum_poly;

}

Polynomial *ply_product(Polynomial *poly1, Polynomial *poly2) {
    Polynomial *p1 = max_deg(poly1, poly2);
    Polynomial *p2 = min_deg(poly1, poly2);
    int deg1 = p1->deg;
    int deg2 = p2->deg;
    int deg = deg1+deg2;

    int k;
    Polynomial *prod_poly = ply_create_poly(deg);
    for (k=0; k<=deg; k++) {
        double sum_k = 0;
            int l;
            for (l=0; l<=k; l++) {
                if ((l<=deg1)&&(k-l<=deg2))
                    sum_k = sum_k+(ply_get_coef(p1, l)*ply_get_coef(p2, k-l));
            }
        ply_set_coef(prod_poly, k, sum_k);
    }

    return prod_poly;
}

Polynomial *ply_scale(double s, Polynomial *p) {
    if (s==0) {
        Polynomial *z = ply_zero();
        return z;
    }
    else if (s!=0) {
        int n = p->deg;
        Polynomial *sp = ply_create_poly(n);

        int i;
        for(i=0; i<=n; i++) {
            double coef = ply_get_coef(p, i);
            ply_set_coef(sp, i, s*coef);
        }

        return sp;
    }
}

// analysis
Polynomial *ply_differentiate(Polynomial *poly, int n) {
    int deg = poly->deg;
    if (deg<n) {
        Polynomial *deriv = ply_create_poly(0);
        ply_set_coef(deriv, 0, 0.0);
        return deriv;
    }
    else {
        int dn_deg = deg-n;
        double dn_coefs[dn_deg+1];
        int i;
        for (i=0; i<=dn_deg; i++) {
            dn_coefs[i] = cnt_factorial(n+i, i+1)*ply_get_coef(poly, n+i);
        }

        while (dn_deg>=0) {
            if (dn_coefs[dn_deg]!=0)
                break;
            if (dn_deg==0)
                break;
            dn_deg = dn_deg-1;
        }

        Polynomial *deriv = ply_create_poly(dn_deg);
        for (i=0; i<=dn_deg; i++) {
            double coef = dn_coefs[i];
            ply_set_coef(deriv, i, coef);
        }

        return deriv;
    }
}

// families of polynomials

// Standard Basis
Polynomial *ply_monomial(int n) {
    Polynomial *p = ply_create_poly(n);
    int i;
    for (i=0; i<=n-1; i++) {
        ply_set_coef(p, i, 0.0);
    }
    ply_set_coef(p, n, 1.0);
    return p;
}

// Legendre (Recursive)
Polynomial *ply_legendre(int n) {
    Polynomial *p0 = ply_monomial(0);
    Polynomial *p1 = ply_monomial(1);

    Polynomial *poly_list[n];
    poly_list[0] = p0;
    poly_list[1] = p1;


    int i;
    for (i=1; i<=n-1; i++) {
        Polynomial *left = ply_scale(2*i+1, ply_product(p1, poly_list[i]));
        Polynomial *right = ply_scale(-i, poly_list[i-1]);
        Polynomial *p = ply_scale(1.0/(i+1), ply_sum(right, left));
        poly_list[i+1] = p;

        ply_delete_poly(left);
        ply_delete_poly(right);
    }

    for(i=0; i<=n-1; i++) {
        ply_delete_poly(poly_list[i]);
    }

    return poly_list[n];
}
