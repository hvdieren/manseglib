#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>

#include <sys/time.h>
#include <time.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

// #define USE_DENSE
#define USE_PRECOND

void iterative_refinement(int n, struct matrix *A, struct matrix *M, DOUBLE *b, DOUBLE *x,
                          int out_maxiter, DOUBLE out_tol, int in_maxiter, DOUBLE in_tol, int step_check,
                          int *out_iter, int *in_iter);

int main(int argc, char *argv[])
{
    if (argc != 7)
    {
        fprintf(stderr, "Missing arguments: algorithm matrix_file out_its out_tol in_its in_tol step_chk\n");
        return 1;
    }

    struct timespec ir_start, ir_end;

    int n, nz;
    struct matrix_coo *coo;
    coo_load(argv[1], &n, &nz, &coo);

#ifdef USE_DENSE
    struct matrix *A = dense_create(n, nz, coo);
#else
    struct matrix *A = csr_create(n, nz, coo);
#endif

    DOUBLE *x = ALLOC(DOUBLE, n); // unknown vector x (what we want to find)
    DOUBLE *b = ALLOC(DOUBLE, n); // well-known vector
    DOUBLE *r = ALLOC(DOUBLE, n); // residual vector
    DOUBLE *s = ALLOC(DOUBLE, n); // randomised starting point (i think)

    vector_rand(n, s);
    // vector_set(n, 1.0 / sqrt(n), s);
    matrix_mult(A, s, b); // b = As

    int out_maxiter = atoi(argv[2]);
    DOUBLE out_tol = atof(argv[3]);
    if (out_tol <= 0.0)
        out_tol = 1.0e-8;

    int in_maxiter = atoi(argv[4]);
    DOUBLE in_tol = atof(argv[5]);
    if (in_tol <= 0.0)
        in_tol = 1.0e-8;
    int step_check = atoi(argv[6]);

    // max error present in matrix
    DOUBLE max = 0;
    for (int i = 0; i < nz; i++)
    {
        DOUBLE x1 = coo[i].a;
        DOUBLE x2 = x1;
        DOUBLE error = (x1 - x2) / x1;
        if (error > max)
            max = error;
    }

#ifdef USE_PRECOND
    struct matrix *M = jacobi_create(n, nz, coo);
    // create jacobian matrix
    // i think this just means we have a diagonal of the matrix coo
    // with all other elements (not on the diagonal) set to zero
#else
    struct matrix *M = NULL;
#endif
    DOUBLE norm = coo_norm_inf(n, nz, coo);
    free(coo);

    printf("# algorithm            : %s\n", argv[0]);
    printf("# sizeof FLOAT         : %lu\n", sizeof(FLOAT));
    // printf("# roundoff: %e\n", epsilon(bits));
    printf("# matrix               : %s\n", argv[1]);
    printf("# problem_size         : %d\n", n);
    printf("# nnz                  : %d\n", nz);
    printf("# matrix_norm          : %e\n", (double)norm);
    printf("# matrix_error         : %e\n", (double)max);
    printf("# bnorm                : %e\n\n", (double)vector_norm2(n, b));
    
    // randomize guesses in x
    vector_rand(n, x);

    int out_iter = 0, in_iter = 0;

    // oprecomp_start();
    // do {
    //

    clock_gettime(CLOCK_MONOTONIC, &ir_start);
    // repeat cg until it converges on a solution or the residual error is too large
    iterative_refinement(n, A, M, b, x, out_maxiter, out_tol, in_maxiter, in_tol, step_check,
                         &out_iter, &in_iter);
    clock_gettime(CLOCK_MONOTONIC, &ir_end);

    //
    // } while (oprecomp_iterate());
    // oprecomp_stop();
    double time_taken = (ir_end.tv_sec - ir_start.tv_sec) * 1e9;
    time_taken = (time_taken + (ir_end.tv_nsec - ir_start.tv_nsec)) * 1e-9;

    matrix_mult(A, x, r); // r = Ax
    vector_xpby(n, b, -1.0, r); // r = x + b^-1
    DOUBLE residual = vector_norm2(n, r);
    DOUBLE normalized_residual = residual / (vector_norm2(n, x) * norm);

    printf("# outer_maxiter        : %d\n", out_maxiter);
    printf("# outer_tolerance      : %e\n", (double)out_tol);
    printf("# inner_maxiter (cg)   : %d\n", in_maxiter);
    printf("# inner_tolerance (cg) : %e\n", (double)in_tol);
    printf("# outer_iterations (ir): %d\n", out_iter);
    printf("# inner_iterations (cg): %d\n", in_iter);
    printf("# residual (ir)        : %e\n", (double)residual);
    printf("# normalized_residual  : %e\n", (double)normalized_residual);

    printf("\n# Time taken           : %.7f s\n", time_taken);

    // printf("\n--- values out ---\n");
    // for(int i = 0; i < n; i++)
    //     printf("%e\n", x[i]);

    return 0;
}