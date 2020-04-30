#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "counting.h"
#include "matrix.h"
#include "polynomial.h"

#define TRUE 1
#define FALSE 0

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
    //PROBLEM corrupted heap
*    double *coef = malloc((deg+1)*sizeof(double)); //PROBLEM
    poly->coefs = coef;

    return poly;
}

void ply_delete_poly(Polynomial *poly) {
    free(poly->coefs);
    free(poly);
}

int ply_get_deg(Polynomial *poly) {
    int deg = poly->deg;
    return deg;
}

double ply_get_coef(Polynomial *poly, int i) {
    return poly->coefs[i];
}

void ply_set_coef(Polynomial *poly, int i, double val) {
    poly->coefs[i] = val;
}

Polynomial *ply_copy(Polynomial *p) {
    Polynomial *copy = ply_create_poly(p->deg);
    int i;
    for (i=0; i<=copy->deg; i++) {
        double coef = ply_get_coef(p, i);
        ply_set_coef(copy, i, coef);
    }

    return copy;
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

// algebra

Polynomial *ply_zero() {
    Polynomial *z = ply_create_poly(0);
    z->coefs[0] = 0.0;
    return z;
    }

int ply_is_zero(Polynomial *p) {
    int deg = p->deg;
    if (ply_get_coef(p, deg)==0)
        return TRUE;
    else
        return FALSE;

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

    //PROBLEM, see ply_create_poly
*    Polynomial *sum_poly = ply_create_poly(deg);//PROBLEM
    for (i=0; i<=deg; i++) {
        double coef = sum_coefs[i];
        ply_set_coef(sum_poly, i, coef);
    }

    return sum_poly;

}

Polynomial *ply_product(Polynomial *poly1, Polynomial *poly2) {
    if (ply_is_zero(poly1)==TRUE)
        return poly1;
    else if (ply_is_zero(poly2)==TRUE)
        return poly2;
    else {
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

static Polynomial **sub_division(Polynomial *f, Polynomial *g) {
    Polynomial *q;
    Polynomial *r;

    int deg_q = f->deg - g->deg;
    q = ply_monomial(deg_q);

    double coef = ply_get_coef(f, f->deg)/ply_get_coef(g, g->deg);
    ply_set_coef(q, deg_q, coef);

    Polynomial *prod = ply_product(g, q);
    Polynomial *scaled = ply_scale(-1.0, prod);
    r = ply_sum(f, scaled);

    ply_delete_poly(prod);
    ply_delete_poly(scaled);

    Polynomial **pair = malloc(2*sizeof(Polynomial *));
    pair[0] = q;
    pair[1] = r;

    return pair;
}

// division algorithm: divides f by g and
// returns struct containing quotient, remainder

// remember to free pair and the polys it contains after use
Polynomial **ply_division(Polynomial *f, Polynomial *g) {
    assert(ply_is_zero(g)==FALSE);
    Polynomial *q0 = ply_zero();
    Polynomial *r0 = ply_copy(f);

    Polynomial *q = q0;
    Polynomial *r = r0;
    Polynomial **pair;

    if (f->deg < g->deg || ply_is_zero(f)==TRUE) {
        pair = malloc(2*sizeof(Polynomial *));
        pair[0] = q;
        pair[1] = r;

        return pair;
    }
    else {
        int change = FALSE;
        while (r->deg >= g->deg && ply_is_zero(r)==FALSE) {
            change = TRUE;
            pair = sub_division(r, g);
            ply_delete_poly(r);
            Polynomial *temp_q = pair[0];
            Polynomial *temp_r = pair[1];

            q = ply_sum(q, temp_q);
            r = ply_copy(temp_r);

            ply_delete_poly(temp_q);
            ply_delete_poly(temp_r);
            free(pair);

        }

        if (change == TRUE)
            ply_delete_poly(q0);

        pair[0] = q;
        pair[1] = r;

        return pair;
    }
}

// return gcd and Bezout coefficients for f,g
Polynomial **ply_gcd(Polynomial *f, Polynomial *g) {
    assert((ply_is_zero(f)==TRUE && ply_is_zero(f)==TRUE)==FALSE);
    Polynomial **result = malloc(3*(sizeof(Polynomial *)));

    if (ply_is_zero(f)==TRUE) {
        result[0] = g;
        result[1] = ply_zero();
        result[2] = ply_monomial(0);

        return result;
    }

    else if (ply_is_zero(g)==TRUE) {
        result[0] = f;
        result[1] = ply_monomial(0);
        result[2] = ply_zero();

        return result;
    }

    else {
        Polynomial *p1 = max_deg(f, g);
        Polynomial *p2 = min_deg(f, g);

        Polynomial **a = malloc(2*sizeof(Polynomial *));
        Polynomial **b = malloc(2*sizeof(Polynomial *));
        Polynomial **c = malloc(2*sizeof(Polynomial *));

        Polynomial *a00 = ply_monomial(0);
        Polynomial *a01 = ply_zero();
        Polynomial *b00 = ply_zero();
        Polynomial *b01 = ply_monomial(0);

        a[0] = a00;
        a[1] = a01;
        b[0] = b00;
        b[1] = b01;

        Polynomial *q0;
        Polynomial *r0;

        Polynomial **qr0 = ply_division(p1, p2);

        q0 = qr0[0];
        r0 = qr0[1];

        Polynomial *q = q0;
        Polynomial *r = r0;

        Polynomial *prod = ply_product(q, b[0]);
        Polynomial *scaled = ply_scale(-1.0, prod);
        //PROBLEM, see ply_sum
*        c[0] = ply_sum(a[0], scaled); //PROBLEM
        ply_delete_poly(prod);
        ply_delete_poly(scaled);

        prod = ply_product(q, b[1]);
        scaled = ply_scale(-1.0, prod);
        c[1] = ply_sum(a[1], scaled);

        ply_delete_poly(prod);
        ply_delete_poly(scaled);

        p1 = p2;
        p2 = r;

        Polynomial **qr;
        int count = 1;
        while (ply_is_zero(r) != TRUE) {
            qr = ply_division(p1, p2);
            q = qr[0];
            r = qr[1];

            ply_delete_poly(a[0]);
            ply_delete_poly(a[1]);
            free(a);

            a = b;
            b = c;

            prod = ply_product(q, b[0]);
            scaled = ply_scale(-1.0, prod);
            c[0] = ply_sum(a[0], scaled);

            ply_delete_poly(prod);
            ply_delete_poly(scaled);

            prod = ply_product(q, b[1]);
            scaled = ply_scale(-1.0, prod);
            c[1] = ply_sum(a[1], scaled);

            ply_delete_poly(prod);
            ply_delete_poly(scaled);

            if (count == 2)
                ply_delete_poly(p1);

            p1 = p2;
            p2 = r;

            if (count < 2)
                count += 1;

            free(qr);
            ply_delete_poly(q);
        }

        double coef = ply_get_coef(p1, p1->deg);
        if (coef != 1.0) {
            Polynomial *t_p1 = p1;
            Polynomial *t_c0 = c[0];
            Polynomial *t_c1 = c[1];

            p1 = ply_scale(1.0/coef, p1);
            c[0] = ply_scale(1.0/coef, c[0]);
            c[1] = ply_scale(1.0/coef, c[1]);

            ply_delete_poly(t_p1);
            ply_delete_poly(t_c0);
            ply_delete_poly(t_c1);
        }

        ply_delete_poly(a00);
        ply_delete_poly(a01);
        ply_delete_poly(b00);
        ply_delete_poly(b01);

        ply_delete_poly(b[0]);
        ply_delete_poly(b[1]);
        free(b);

        ply_delete_poly(q0);
        ply_delete_poly(r0);

        free(qr0);
        ply_delete_poly(r);

        result[0] = p1;
        result[1] = c[0];
        result[2] = c[1];

        return result;

    }
}

//---------------------------------------------------------------------
/*

Polynomial *ply_modulo(Polynomial *p, Polynomial *modulus) {

}
*/





// analysis

// Horner's method
double ply_evaluate(Polynomial *poly, double x) {
    int deg = poly->deg;
    double coefn = ply_get_coef(poly, deg);
    assert(coefn!=0);

    double val = coefn;
    if (deg==0) {
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
