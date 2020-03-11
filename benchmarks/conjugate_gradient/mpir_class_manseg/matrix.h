#ifndef __cg_MATRIX_H__
#define __cg_MATRIX_H__

// #include "../../../mantissaSegmentation.hpp"

// definitions for matrices
class matrix
{
public:
    int n;
    
    void (*dmult)(matrix *, PairsArray *, PairsArray *);
    void (*smult)(matrix *, HeadsArray *, HeadsArray *);
};

class matrix_coo : public matrix
{
public:
    int *i, *j;
    // double *a;
    ManSegArray *a;
};

class matrix_csr : public matrix
{
public:
    int n;
    int *i, *j;
    // double *A;
    ManSegArray *A;
};

class matrix_dense : public matrix
{
public:
    int n;
    // double *A;
    ManSegArray *A;
};

class precond_jacobi : public matrix
{
public:
    int n;
    // float *d;
    ManSegArray *d;
};

static inline void matrix_mult(matrix *mat, PairsArray *x, PairsArray *y) {
    mat->dmult(mat, x, y);
}

static inline void floatm_mult(matrix *mat, HeadsArray *x, HeadsArray *y) {
    mat->smult(mat, x, y);
}

extern matrix_coo *coo_load(const char *fname, int *n, int *nz);
extern double coo_norm_inf(int n, int nz, matrix_coo *coo);
extern double coo_max_nz(int n, int nz, matrix_coo *coo);

extern matrix *csr_create(int n, int nz, matrix_coo *coo);
extern matrix *dense_create(int n, int nz, matrix_coo *coo);
extern matrix *jacobi_create(int n, int nz, matrix_coo *coo);


#endif
