#ifndef NUMBER_FIELD_H
#define NUMBER_FIELD_H

#include "matrix.h"

// number field
typedef struct AlgebraicNumber AlgebraicNumber;

struct AlgebraicNumber {
    Polynomial min_poly;
    Polynomial number;
};

// structure
AlgebraicNumber nf_create(Polynomial min_poly, Polynomial number);
void nf_print(AlgebraicNumber x);

// algebra
AlgebraicNumber nf_sum(AlgebraicNumber x, AlgebraicNumber y);
AlgebraicNumber nf_neg(AlgebraicNumber x);
AlgebraicNumber nf_product(AlgebraicNumber x, AlgebraicNumber y);
AlgebraicNumber nf_inv(AlgebraicNumber x);



#endif
