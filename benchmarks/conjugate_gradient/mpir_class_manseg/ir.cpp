#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter,
						FLOAT *x_prev, FLOAT2 *component_diff_x_cur, FLOAT2 *component_diff_x_prev);

void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter);

void iterative_refinement(int n, matrix *A, matrix *M, DOUBLE *b, DOUBLE *x, int out_maxiter, DOUBLE out_tol, 
    int in_maxiter, DOUBLE in_tol, int step_check, int *out_iter, int *in_iter)
{
    DOUBLE *e = ALLOC(DOUBLE, n);
    FLOAT *r = ALLOC(FLOAT, n);
    FLOAT *d = ALLOC(FLOAT, n);
	// FLOAT *d_prev = ALLOC(FLOAT, n);
	// DOUBLE *component_diff_x_cur = ALLOC(DOUBLE, n);
	// DOUBLE *component_diff_x_prev = ALLOC(DOUBLE, n);

    mixed_copy(n, x, d); // d = (float)x

	// floatm_copy(n, d, d_prev);

    vector_set(n, 0.0, x);

	// vector_set(n, 0.0f, component_diff_x_cur);
	// vector_set(n, 0.0f, component_diff_x_prev);

    mixed_copy(n, b, r); // r = (float)b

    DOUBLE residual;
    do
    {
        // conjugate_gradient(n, A, M, r, d, in_maxiter, in_tol, step_check, in_iter,
		// 					d_prev, component_diff_x_cur, component_diff_x_prev);
		conjugate_gradient(n, A, M, r, d, in_maxiter, in_tol, step_check, in_iter);
        
        mixed_axpy(n, 1.0, d, x); // x = x + d
        matrix_mult(A, x, e);
        vector_xpby(n, b, -1.0, e); // r = b - Ax
        residual = vector_norm2(n, e);

        mixed_copy(n, e, r);
        floatm_set(n, 0.0, d);
		// floatm_set(n, 0.0f, d_prev);

        // *energy += iter * (bits + 12) / 8;
        // printf("%d %d %d %e %e\n", *in_iter, 0, bits, (double)residual, (double)residual);
        (*out_iter)++;		
    } while ((residual > out_tol) && (*out_iter < out_maxiter));

	// *in_iter += *out_iter;

    if (residual <= out_tol)
        printf("\n==== residual less than out_tol ====\n");

    FREE(e);
    FREE(r);
    FREE(d);
	// FREE(d_prev);
	// FREE(component_diff_x_cur);    FREE(component_diff_x_prev);
}