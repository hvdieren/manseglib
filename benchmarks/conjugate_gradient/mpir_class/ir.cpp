#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter);

void iterative_refinement(int n, matrix *A, matrix *M, DOUBLE *b, DOUBLE *x, int out_maxiter, DOUBLE out_tol, 
    int in_maxiter, DOUBLE in_tol, int step_check, int *out_iter, int *in_iter)
{
    DOUBLE *e = ALLOC(DOUBLE, n);
    FLOAT *r = ALLOC(FLOAT, n);
    FLOAT *d = ALLOC(FLOAT, n);

    mixed_copy(n, x, d); // d = (float)x
    vector_set(n, 0.0, x);
    mixed_copy(n, b, r); // r = (float)b

    DOUBLE residual;
    do
    {
		conjugate_gradient(n, A, M, r, d, in_maxiter, in_tol, step_check, in_iter);
        
        mixed_axpy(n, 1.0, d, x); // x = x + d
		matrix_mult(A, x, e);
        vector_xpby(n, b, -1.0, e); // r = b - Ax

        residual = vector_norm2(n, e);
        mixed_copy(n, e, r);
        floatm_set(n, 0.0, d);

		printf("%d: outer: residual = %e\n", *out_iter, residual);

        // *energy += iter * (bits + 12) / 8;
        // printf("%d %d %d %e %e\n", *in_iter, 0, bits, (double)residual, (double)residual);

        (*out_iter)++;		
    } while ((residual > out_tol) && (*out_iter < out_maxiter));

    if (residual <= out_tol)
        printf("\n==== residual less than out_tol ====\n");
	if(*out_iter >= out_maxiter)
		printf("out_iter >= out_maxiter\n");

    FREE(e);
    FREE(r);
    FREE(d);
}