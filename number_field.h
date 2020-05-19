#ifndef NUMBER_FIELD_H
#define NUMBER_FIELD_H

#include "matrix.h"

// number field
typedef struct NFNumber NFNumber;

struct NFNumber {
    Polynomial min_poly;
    Polynomial number;
};

// structure
NFNumber nf_create(Polynomial min_poly, Polynomial number);
void nf_print(NFNumber x);

// algebra
NFNumber nf_sum(NFNumber x, NFNumber y);
NFNumber nf_neg(NFNumber x);
NFNumber nf_product(NFNumber x, NFNumber y);
NFNumber nf_inv(NFNumber x);



#endif
