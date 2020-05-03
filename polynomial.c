#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "counting.h"
#include "polynomial.h"
#include "matrix.h"

#define TRUE 1
#define FALSE 0

// helpful functions
static Polynomial max_deg(Polynomial p1, Polynomial p2) {
    int deg1 = p1.deg;
    int deg2 = p2.deg;
    if (deg1>deg2)
        return p1;
    else if (deg2>deg1)
        return p2;
    else {
        if (ply_is_zero(p1))
            return p2;
        else
            return p1;
    }
}

static Polynomial min_deg(Polynomial p1, Polynomial p2) {
    int deg1 = p1.deg;
    int deg2 = p2.deg;
    if (deg1>deg2)
        return p2;
    else if (deg2>deg1)
        return p1;
    else {
        if (ply_is_zero(p1))
            return p1;
        else
            return p2;
    }
}

// data structure
Polynomial ply_create(int deg) {
    assert(deg >= 0 && deg <= MAX_DEG);
    Polynomial poly;
    poly.deg = deg;

    return poly;
}

int ply_get_deg(Polynomial poly) {
    int deg = poly.deg;
    return deg;
}

double ply_get_coef(Polynomial poly, int i) {
    return poly.coefs[i];
}

void ply_set_coef(Polynomial *poly, int i, double val) {
    poly->coefs[i] = val;
}

Polynomial ply_copy(Polynomial p) {
    Polynomial copy = ply_create(p.deg);
    int i;
    for (i=0; i<=copy.deg; i++) {
        double coef = ply_get_coef(p, i);
        ply_set_coef(&copy, i, coef);
    }

    return copy;
}

void ply_print(Polynomial poly) {
    int deg = poly.deg;
    double coef0 = ply_get_coef(poly, 0);
    double coefn = ply_get_coef(poly, deg); //PROBLEM
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

int ply_equal(Polynomial p, Polynomial q) {
    int m = p.deg;
    int n = q.deg;

    if (m != n)
        return FALSE;
    else {
        int i;
        for (i=0; i<=m; i++){
            double co_p = ply_get_coef(p, i);
            double co_q = ply_get_coef(q, i);
            if (co_p != co_q)
                return FALSE;
        }

        return TRUE;
    }
}

Polynomial ply_zero() {
    Polynomial z = ply_create(0);
    z.coefs[0] = 0.0;
    return z;
    }

int ply_is_zero(Polynomial p) {
    int deg = p.deg;
    if (ply_get_coef(p, deg)==0)
        return TRUE;
    else
        return FALSE;

}

Polynomial ply_sum(Polynomial poly1, Polynomial poly2) {

    Polynomial p1 = max_deg(poly1, poly2);
    Polynomial p2 = min_deg(poly1, poly2);
    int deg1 = p1.deg;
    int deg2 = p2.deg;

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

    Polynomial sum_poly = ply_create(deg);
    for (i=0; i<=deg; i++) {
        double coef = sum_coefs[i];
        ply_set_coef(&sum_poly, i, coef);
    }

    return sum_poly;

}

Polynomial ply_product(Polynomial poly1, Polynomial poly2) {
    if (ply_is_zero(poly1)==TRUE)
        return poly1;
    else if (ply_is_zero(poly2)==TRUE)
        return poly2;
    else {
        Polynomial p1 = max_deg(poly1, poly2);
        Polynomial p2 = min_deg(poly1, poly2);
        int deg1 = p1.deg;
        int deg2 = p2.deg;
        int deg = deg1+deg2;

        int k;
        Polynomial prod_poly = ply_create(deg);
        for (k=0; k<=deg; k++) {
            double sum_k = 0;
                int l;
                for (l=0; l<=k; l++) {
                    if ((l<=deg1)&&(k-l<=deg2))
                        sum_k = sum_k+(ply_get_coef(p1, l)*ply_get_coef(p2, k-l));
                }
            ply_set_coef(&prod_poly, k, sum_k);
        }

        return prod_poly;
    }
}

Polynomial ply_scale(double s, Polynomial p) {
    if (s==0) {
        Polynomial z = ply_zero();
        return z;
    }
    else {
        int n = p.deg;
        Polynomial sp = ply_create(n);

        int i;
        for(i=0; i<=n; i++) {
            double coef = ply_get_coef(p, i);
            ply_set_coef(&sp, i, s*coef);
        }

        return sp;
    }
}











//---------------------------------------------------------------
static PolyMatrix sub_division(Polynomial f, Polynomial g) {
    Polynomial q;
    Polynomial r;

    int deg_q = f.deg - g.deg;
    q = ply_monomial(deg_q);

    double coef = ply_get_coef(f, f.deg)/ply_get_coef(g, g.deg);
    ply_set_coef(&q, deg_q, coef);

    Polynomial prod = ply_product(g, q);
    Polynomial scaled = ply_scale(-1.0, prod);
    r = ply_sum(f, scaled);

    PolyMatrix pair = pymat_create(1, 2);
    pymat_set_element(&pair, 0, 0, q);
    pymat_set_element(&pair, 0, 1, r);

    return pair;
}

// division algorithm: divides f by g and
// returns struct containing quotient, remainder

// remember to free pair and the polys it contains after use
PolyMatrix ply_division(Polynomial f, Polynomial g) {
    printf("line 277 \n" );
    assert(ply_is_zero(g)==FALSE);
    Polynomial q0 = ply_zero();
    Polynomial r0 = ply_copy(f);

    Polynomial q = q0;
    Polynomial r = r0;
    PolyMatrix pair;

    printf("line 285\n" );

    if (f.deg < g.deg || ply_is_zero(f)==TRUE) {
        pair = pymat_create(1,2);
        pymat_set_element(&pair, 0, 0, q);
        pymat_set_element(&pair, 0, 1, r);

        return pair;
    }
    else {
        while (r.deg >= g.deg && ply_is_zero(r)==FALSE) {
            pair = sub_division(r, g);
            Polynomial temp_q = pymat_get_element(pair, 0, 0);
            Polynomial temp_r = pymat_get_element(pair, 0, 1);

            q = ply_sum(q, temp_q);
            r = ply_copy(temp_r);

        }

        pair = pymat_create(1,2);
        pymat_set_element(&pair, 0, 0, q);
        pymat_set_element(&pair, 0, 1, r);

        return pair;
    }
}

// return gcd and Bezout coefficients for f,g
PolyMatrix ply_gcd(Polynomial f, Polynomial g) {
    assert((ply_is_zero(f)==TRUE && ply_is_zero(f)==TRUE)==FALSE);
    PolyMatrix result = pymat_create(1, 3);

    if (ply_is_zero(f)==TRUE) {
        pymat_set_element(&result, 0, 0, g);
        pymat_set_element(&result, 0, 1, ply_zero());
        pymat_set_element(&result, 0, 2, ply_monomial(0));

        return result;
    }

    else if (ply_is_zero(g)==TRUE) {
        pymat_set_element(&result, 0, 0, f);
        pymat_set_element(&result, 0, 1, ply_monomial(0));
        pymat_set_element(&result, 0, 2, ply_zero());

        return result;
    }

    else {
        printf("line 333\n" );
        Polynomial p1 = max_deg(f, g);
        Polynomial p2 = min_deg(f, g);

        PolyMatrix a = pymat_create(1, 2);
        PolyMatrix b = pymat_create(1, 2);
        PolyMatrix c;

        Polynomial a00 = ply_monomial(0);
        Polynomial a01 = ply_zero();
        Polynomial b00 = ply_zero();
        Polynomial b01 = ply_monomial(0);

        pymat_set_element(&a, 0, 0, a00);
        pymat_set_element(&a, 0, 1, a01);
        pymat_set_element(&b, 0, 0, b00);
        pymat_set_element(&b, 0, 1, b01);

        while(1) {
            printf("line 353\n" );
            PolyMatrix qr = ply_division(p1, p2); \\problem
            printf("line 354\n" );
            Polynomial q = pymat_get_element(qr, 0, 0);
            Polynomial r = pymat_get_element(qr, 0, 1);

            if (ply_is_zero(r))
                break;

            PolyMatrix scaled = pymat_poly_scale(q, b);
            PolyMatrix neg_qb = pymat_scale(-1.0, scaled);
            c = pymat_sum(a, neg_qb);

            p1 = ply_copy(p2);
            p2 = ply_copy(r);

            a = b;
            b = c;

        }
        printf("line 373\n" );
        double coef = ply_get_coef(p2, p2.deg);
        if (coef != 1.0) {
            Polynomial t_p2 = p2;
            PolyMatrix t_c = c;

            p2 = ply_scale(1.0/coef, t_p2);
            c = pymat_scale(1.0/coef, t_c);

        }

        Polynomial c00 = pymat_get_element(c, 0, 0);
        Polynomial c01 = pymat_get_element(c, 0, 1);
        pymat_set_element(&result, 0, 0 ,p2);
        pymat_set_element(&result, 0, 1, c00);
        pymat_set_element(&result, 0, 1, c01);

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
double ply_evaluate(Polynomial poly, double x) {
    int deg = poly.deg;
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

Polynomial ply_differentiate(Polynomial poly, int n) {
    int deg = poly.deg;
    if (deg<n) {
        Polynomial deriv = ply_create(0);
        ply_set_coef(&deriv, 0, 0.0);
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

        Polynomial deriv = ply_create(dn_deg);
        for (i=0; i<=dn_deg; i++) {
            double coef = dn_coefs[i];
            ply_set_coef(&deriv, i, coef);
        }

        return deriv;
    }
}

// families of polynomials

// Standard Basis
Polynomial ply_monomial(int n) {
    Polynomial p = ply_create(n);
    int i;
    for (i=0; i<=n-1; i++) {
        ply_set_coef(&p, i, 0.0);
    }
    ply_set_coef(&p, n, 1.0);
    return p;
}

// Legendre (Recursive)
Polynomial ply_legendre(int n) {
    Polynomial p0 = ply_monomial(0);
    Polynomial p1 = ply_monomial(1);

    Polynomial poly_list[n];
    poly_list[0] = p0;
    poly_list[1] = p1;

    int i;
    for (i=1; i<=n-1; i++) {

        Polynomial prod = ply_product(p1, poly_list[i]);
        Polynomial left = ply_scale(2*i+1, prod);


        Polynomial right = ply_scale(-i, poly_list[i-1]);
        Polynomial sum = ply_sum(right, left);
        Polynomial p = ply_scale(1.0/(i+1), sum);

        poly_list[i+1] = p;
    }

    return poly_list[n];
}
