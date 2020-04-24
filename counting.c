// counting

#include <stdio.h>
#include <assert.h>
#include <inttypes.h>

int64_t cnt_factorial(int64_t n, int64_t k) {
    assert(k<=n);
    assert(0<=k);
    if (n==0) {
        return 1;
    }
    else {
        assert(k>=1);
        int64_t prod = 1;
        int i;
        for (i=k; i<=n; i++) {
            prod = prod*i;
        }
        return prod;
    }
}

int64_t cnt_permutation(int n, int k) {
    assert(0<=k);
    if (k==0) {
        return 1;
    }
    if (n<k) {
        return 0;
    }
    else {
        return cnt_factorial(n,n-k+1);
    }
}

int64_t cnt_combination(int n, int k) {
    assert(0<=k);
    if (k==0) {
        return 1;
    }
    if (n<k) {
        return 0;
    }
    else {
        return cnt_permutation(n,k)/cnt_factorial(k,1);
    }
}
