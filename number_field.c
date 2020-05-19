#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "polynomial.h"
#include "number_field.h"

#define TRUE 1
#define FALSE 0

// structure

NFNumber nf_create(Polynomial min_poly, Polynomial number) {
    NFNumber x;
    x.min_poly = min_poly;
    if (number.deg > min_poly.deg) {
        x.number = pymod_reduce(number, min_poly);
    }
    else {
        x.number = number;
    }


    return x;
}

void nf_print(NFNumber x) {
    Polynomial poly = x.number;
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
            printf("%.2fa^%d ", coef, k);
            for (i=k+1; i<=deg-1; i++) {
                double coef = ply_get_coef(poly, i);
                if (coef>0) {
                    printf("+ %.2fa^%d ", coef, i);
                }
                else if (coef<0)
                    printf("- %.2fa^%d ", -coef, i);
                }
            if (coefn>0)
                printf("+ %.2fa^%d\n", coefn, deg);
            else
                printf("- %.2fa^%d\n", -coefn, deg);
        }

        else if (k==deg) {
            printf("%.3fa^%d\n", coefn, deg);
        }
    }

    else if (deg==0) {
        printf("%.2f\n", coef0);
    }
}
//-----------------------------------------------

// algebra

NFNumber nf_sum(NFNumber x, NFNumber y) {
    Polynomial z_num = pymod_sum(x.number, y.number, x.min_poly);
    NFNumber z = nf_create(x.min_poly, z_num);

    return z;
}

NFNumber nf_neg(NFNumber x) {
    Polynomial min_poly = x.min_poly;
    Polynomial neg_num = ply_neg(x.number);

    NFNumber neg_x = nf_create(min_poly, neg_num);

    return neg_x;
}

NFNumber nf_product(NFNumber x, NFNumber y) {
    Polynomial z_num = pymod_product(x.number, y.number, x.min_poly);
    NFNumber z = nf_create(x.min_poly, z_num);

    return z;
}

NFNumber nf_inv(NFNumber x) {
    assert(ply_is_zero(x.number)==FALSE);
    Polynomial inv_num = pymod_inv(x.number, x.min_poly);

    NFNumber inv = nf_create(x.min_poly, inv_num);

    return inv;
}
