#include <stdlib.h>

// helpful functions
static int min(double a, double b) {
    if (a < b)
        return a;
    else
        return b;
}
//----------------------------------------------------------------------------

// array algorithms

void arr_dubswap(double *arr, int i, int j) {
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

void arr_intswap(int *arr, int i, int j) {
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

double arr_sum(double *arr, int length) {
    double sum = 0;
    int i;
    for (i=0; i<length; i=i+1)
        sum = sum + arr[i];
    return sum;
}

// merge-sort an array

// merge two sorted subarrays arr[l,l+1,...m] arr[m+1,m+2,...,r]
static void merge(double *arr, int l, int m, int r) {
    int i, j, k;
    int n1 = m-l+1;
    int n2 = r-m;

    //create temporary left and right arrays
    double la[n1], ra[n2];

    // copy data to temp arrays
    for (i=0; i<n1; i=i+1)
        la[i] = arr[l+i];

    for (j=0; j<n2; j=j+1)
        ra[j] = arr[m+1+j];

    // merge the temp arrays
    i = 0;
    j = 0;
    k = l;

    while (i<n1 && j<n2) {
        if (la[i] <= ra[j]){
            arr[k] = la[i];
            i = i+1;
        }
        else {
            arr[k] = ra[j];
            j = j+1;
        }
        k = k+1;
    }

    while (i<n1) {
        arr[k] = la[i];
        i = i+1;
        k = k+1;
    }

    while (j<n2) {
        arr[k] = ra[j];
        j = j+1;
        k = k+1;
    }
}

void arr_merge_sort(double *arr, int l, int r) {
    if (l < r) {
        int m = l + (r-l)/2;
        arr_merge_sort(arr, l, m);
        arr_merge_sort(arr, m+1, r);

        merge(arr, l, m, r);
    }
}

//quick-sort an array
typedef struct Tuple_ Tuple;

struct Tuple_ {
    int values[2];
};

// parttion arr with pivot as last element
static Tuple partition(double *arr, int l, int r) {
    int i = l, j = l, k = 0;
    while (j < r-k) {
        if (arr[j] == arr[r-k]) {
            if (arr[j] == arr[r-k-1]) {
                k = k+1;
                continue;
            }
            else {
                arr_dubswap(arr, j, r-k-1);
                k = k+1;
            }
        }
        if (arr[j] < arr[r-k]) {
            if (i != j)
                arr_dubswap(arr, i, j);
            i = i+1;
        }
        j = j+1;
    }

    int n, m = min(j-i-1, k);
    for (n=0; n<=m; n=n+1)
        arr_dubswap(arr, i+n, r-n);

        Tuple x;
    if (k < r-l) {
        x.values[0] = i-1;
        x.values[1] = i+k+1;
    }
    else {
        x.values[0] = l;
        x.values[1] = r;
    }

    return x;
}

void arr_quick_sort(double *arr, int l, int r) {
    if (l < r) {
        Tuple x = partition(arr, l, r);
        arr_quick_sort(arr, l, x.values[0]);
        arr_quick_sort(arr, x.values[1], r);
    }
}

// binary-search an array

//search

int arr_binary_search(double *arr, int l, int r, int x) {
    if (r > l) {
        int m = l + (r-l)/2;

        if (arr[m] == x)
            return m;
        else if (arr[m]<x)
            return arr_binary_search(arr, m+1, r, x);
        else
            return arr_binary_search(arr, l, m, x);
    }
    else {
        if (arr[l] == x)
            return l;
        else
            return -1;
    }
}
