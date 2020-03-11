
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

#include "mmio.h"
#include <tgmath.h> // interferes with mmio.h

#include "cg.h"
#include "matrix.h"

int partition_coo(matrix_coo *coo, int lo, int hi)
{
    int pivot = coo->i[(hi + lo)/2];
    int i = lo - 1;
    int j = hi + 1;

    while(true)
    {
        while(coo->i[++i] < pivot);
        while(coo->i[--j] > pivot);
        if(i >= j)
            return j;
        std::swap(coo->i[i], coo->i[j]);
        std::swap(coo->j[i], coo->j[j]);

        // std::swap(coo->a[i], coo->a[j]);
        double temp = coo->a->pairs[j];
        coo->a->pairs[j] = coo->a->pairs[i];
        coo->a->pairs[i] = temp;
    }
}

void quicksort_coo(matrix_coo *coo, int lo, int hi)
{
    if(lo < hi)
    {
        int p = partition_coo(coo, lo, hi);
        quicksort_coo(coo, lo, p);
        quicksort_coo(coo, (p + 1), hi);
    }
}

matrix_coo* coo_load(const char *fname, int *n, int *nz)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Error opening file: %s\n", fname);
        exit(1);
    }
    MM_typecode matcode;
    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner\n");
        exit(1);
    }
    if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) || mm_is_complex(matcode) || !mm_is_symmetric(matcode)) {
        fprintf(stderr, "This application does not support the Market Market type: %s\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }
    int M, N, NZ;
    if (mm_read_mtx_crd_size(f, &M, &N, &NZ) != 0) {
        fprintf(stderr, "Could not parse matrix size\n");
        exit(1);
    }
    if (M != N) {
        fprintf(stderr, "Matrix is not square\n");
        exit(1);
    }

    matrix_coo* coo = new matrix_coo();
    coo->i = new int[2*NZ] ();
    coo->j = new int[2*NZ] ();
    // coo->a = new double[2*NZ] ();
    coo->a = new ManSegArray(2*NZ);

    int k = 0;
    for (int l = 0; l < NZ; l++) {
        double real;
        if (mm_read_mtx_crd_entry(f, &coo->i[k], &coo->j[k], &real, NULL, matcode) != 0) {
            fprintf(stderr, "Error reading matrix element %i\n", l);
            exit(1);
        }
        coo->i[k]--;
        coo->j[k]--;
        coo->a->pairs[k] = real;
        if (coo->i[k] == coo->j[k]) k++;
        else {
            coo->i[k + 1] = coo->j[k];
            coo->j[k + 1] = coo->i[k];
            coo->a->pairs[k + 1] = coo->a->pairs[k];
            k += 2;
        }
    }
    fclose(f);

    *n = N; *nz = k;
    quicksort_coo(coo, 0, k-1);

    return coo;
}

double coo_norm_inf(int n, int nz, matrix_coo *coo)
{
    double *amax = CALLOC(double, n);
    for (int i = 0; i < nz; i++) {
        int row = coo->i[i];
        int val = coo->a->pairs[i];
        amax[row] += abs(val);
    }
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        if (amax[i] > norm) norm = amax[i];
    }
    FREE(amax);
    return norm;
}

double coo_max_nz(int n, int nz, matrix_coo *coo)
{
    double *m = CALLOC(double, n);
    for (int i = 0; i < nz; i++) {
        int row = coo->i[i];
        m[row]++;
    }
    int r = 0.0;
    for (int i = 0; i < n; i++) {
        if (m[i] > r) r = m[i];
    }
    FREE(m);
    return r;
}

// CSR matrix
void csr_dmult(matrix_csr *mat, PairsArray *x, PairsArray *y)
{
    for (int k = 0; k < mat->n; k++) {
        DOUBLE t = 0.0;
        for (int l = mat->i[k]; l < mat->i[k + 1]; l++)
            t += mat->A->pairs[l] * (*x)[mat->j[l]];
        (*y)[k] = t;
    }
}

void csr_smult(matrix_csr *mat, HeadsArray *x, HeadsArray *y)
{
    for (int k = 0; k < mat->n; k++) {
        FLOAT2 t = 0.0;
        for (int l = mat->i[k]; l < mat->i[k + 1]; l++)
            t += mat->A->heads[l] * (*x)[mat->j[l]];
        (*y)[k] = t;
    }
}

matrix *csr_create(int n, int nz, matrix_coo *coo)
{
    int *i = ALLOC(int, n + 1);
    int *j = ALLOC(int, nz);
    // DOUBLE *A = ALLOC(DOUBLE, nz);
    ManSegArray *A = new ManSegArray(nz);

    i[0] = 0;
    int l = 0;
    for (int k = 0; k < n; k++) {
        while (l < nz && coo->i[l] == k) {
            j[l] = coo->j[l];
            A->pairs[l] = coo->a->pairs[l];
            l++;
        }
        i[k + 1] = l;
    }

    matrix_csr *mat = new matrix_csr();
    mat->n = n;
    mat->i = i;
    mat->j = j;
    mat->A = A;
    mat->dmult = (void (*)(matrix *, PairsArray *, PairsArray *))csr_dmult;
    mat->smult = (void (*)(matrix *, HeadsArray *, HeadsArray *))csr_smult;
    return (matrix *)mat;
}

// dense matrix
void dense_dmult(matrix_dense *mat, PairsArray *x, PairsArray *y)
{
    for (int i = 0; i < mat->n; i++) {
        DOUBLE t = 0.0;
        for (int j = 0; j < mat->n; j++)
            t += mat->A->pairs[i * mat->n + j] * (*x)[j];
        (*y)[i] = t;
    }
}

void dense_smult(uint8_t m, matrix_dense *mat, HeadsArray *x, HeadsArray *y)
{
    for (int i = 0; i < mat->n; i++) {
        FLOAT2 t = 0.0;
        for (int j = 0; j < mat->n; j++)
            t += mat->A->heads[i * mat->n + j] * (*x)[j];
        (*y)[i] = t;
    }
}

matrix *dense_create(int n, int nz, matrix_coo *coo)
{
    // DOUBLE *A = CALLOC(DOUBLE, n * n);
    ManSegArray *A = new ManSegArray(n*n);

    // row major format
    for (int k = 0; k < nz; k++) A->pairs[coo->i[k] * n + coo->j[k]] = coo->a->pairs[k];

    matrix_dense *mat = new matrix_dense();
    mat->n = n;
    mat->dmult = (void (*)(matrix *, PairsArray *, PairsArray *))dense_dmult;
    mat->smult = (void (*)(matrix *, HeadsArray *, HeadsArray *))dense_smult;
    mat->A = A;
    return (matrix *)mat;
}

// Jacobi preconditioner
void jacobi_smult(precond_jacobi *pre, HeadsArray *x, HeadsArray *y)
{
    for (int k = 0; k < pre->n; k++) {
       (*y)[k] = (*x)[k] / pre->d->heads[k];
    }
}

struct matrix *jacobi_create(int n, int nz, matrix_coo *coo)
{
    // FLOAT *d = ALLOC(FLOAT, n);
    ManSegArray *d = new ManSegArray(n);
    for (int k = 0; k < n; k++) d[k] = 0.0;
    for (int k = 0; k < nz; k++)
        if (coo->i[k] == coo->j[k])
            d->heads[coo->i[k]] = coo->a->heads[k];

    precond_jacobi *pre = new precond_jacobi();
    pre->n = n;
    pre->dmult = NULL;
    pre->smult = (void (*)(matrix *, HeadsArray *, HeadsArray *))jacobi_smult;
    pre->d = d;
    return (matrix *)pre;
}