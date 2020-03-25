#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

void inline check_x_delta(int *n, matrix *A, FLOAT *x, FLOAT *x_prev, int *in_iter, FLOAT2 *residual)
{
	// DOUBLE x_diff = floatm_max_diff(*n, x_prev, x);							// find max delta in solution (x) vector
	// printf("x_diff=%e\n", x_diff);
	// if (x_diff <= AdaptivePrecisionBound) {									// if close to the max precision of the reduced format, switch to use full precision
	// 	mat_increase_precision(A);
	// 	printf("increased precision at iter %d, x_diff=%e, residual=%e\n", *in_iter, x_diff, *residual);
	// }
}

// void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter,
// 						FLOAT *x_prev, FLOAT2 *component_diff_x_cur, FLOAT2 *component_diff_x_prev)
void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter)
{
    int iter = 0;
    FLOAT2 alpha, beta;
	FLOAT2 rho; // rho = r(k-1)^Tz(k-1)
	FLOAT2 tau; // tau = r(k)^Tz(k)
	FLOAT2 tol; // tol = ||r(k)||v2
	const DOUBLE eps = pow(2.0, -52.0);
	const int check_freq = step_check/10;
	const double switchgrad = 2;

	// x(0) = [0, 0, .., 0]

    FLOAT *r = ALLOC(FLOAT, n);
    FLOAT *p = ALLOC(FLOAT, n); // this may be the actual "conjugate gradient" (conjugate vectors)
    FLOAT *z = ALLOC(FLOAT, n);
	// FLOAT *w = ALLOC(FLOAT, n);

    FLOAT *tr = ALLOC(FLOAT, n);

	/* FLOAT *x_prev = ALLOC(FLOAT, n);
	floatm_copy(n, x, x_prev);

	FLOAT2 *component_diff_x_cur = CALLOC(FLOAT2, n);
	FLOAT2 *component_diff_x_prev = CALLOC(FLOAT2, n); */

	/* FLOAT *c = ALLOC(FLOAT, n);
	FLOAT *cphi = ALLOC(FLOAT, n);

	FLOAT *t1 = ALLOC(FLOAT, n);
	FLOAT *t2 = ALLOC(FLOAT, n);

	int *tcond = ALLOC(int, n); */

	// r(0) = b - Ax(0)
    floatm_mult(A, x, r);
    floatm_xpby(n, b, -1.0, r);

    if (M) {
		// Mz(0) = r(0)
        floatm_mult(M, r, z);
		// calculate rho
        rho = floatm_dot(n, r, z);
		// calculate tol
		tol = floatm_norm2(n, r);
    } else { // small optimization
		// z(0) = r(0)
        floatm_copy(n, r, z);
		// rho = r(0)^Tr(0)
        rho = floatm_dot(n, r, r);
        // tol = sqrt(rho)
		tol = sqrt(rho);
    }
	// p(1) = z(0)
	floatm_copy(n, z, p);

    int step = 0;
    FLOAT2 residual = tol;

    while ((iter < maxiter) && (tol > umbral)) {
		// w = Ap(k+1)
		floatm_mult(A, p, z);

		if (step < step_check) {
			step++;
		}
		else {
			floatm_mult(A, x, tr);		 // tr = Ax
			floatm_xpby(n, b, -1.0, tr); // tr = b - Ax
			residual = floatm_norm2(n, tr);
			printf("# rescheck: total_cg_iter=%d current_iter=%d tol=%e resid=%e\n", *in_iter, iter, (double)tol, (double)residual);

			double res_over_tol = residual/tol;

			if(!A->useTail && (res_over_tol > switchgrad)) 
				mat_increase_precision(A);

			printf("%d: %e\n", *in_iter, residual);

			if (res_over_tol > 10) {
				break;
			}
			step = 1;
		}

		// alpha(k+1) = rho / p(k+1)w
		alpha = rho / floatm_dot(n, z, p);

		// x(k+1) = x + alpha(k+1)p(k+1)
		floatm_axpy(n, alpha, p, x);

		// r(k+1) = r(k) - alpha(k+1)(w)
		floatm_axpy(n, -alpha, z, r);

		// printf("|r(i)|^2 = [");
		// for(int i = 0; i < n; i++)
		// 	printf("%e ", r[i]*r[i]);
		// printf("]\n");

		// apply preconditioner
		if (M) {
			// Solve Mz(k) = r(k)
		    floatm_mult(M, r, z);
			// calculate tau
		    tau = floatm_dot(n, r, z);
			// calculate tol
		    tol = floatm_norm2(n, r);
		} else {
			// tau = inner product(r)
		    tau = floatm_dot(n, r, r);
			// tol = sqrt(tau)
		    tol = sqrt(tau);
		}
		
		// beta = tau / rho
		beta =  tau / rho;
		// update r(k-1)^Tz(k-1) to r(k)^Tz(k)
		rho = tau;
		
		// p(k+1) = z(k) + beta(k) + p(k)
		if (M) {
			floatm_xpby(n, z, beta, p);
		}
		else {
			floatm_xpby(n, r, beta, p);
		}

		// if((iter % check_freq) == 0) {
		// 	floatm_max_abs_diff(n, x, x_prev, component_diff_x_cur, eps);
		// 	// floatm_ratio(n, x, x_prev, component_ratio_x);

		// 	FLOAT2 min_abs_diff = __DBL_MAX__;
		// 	for(int i = 0; i < n; i++)
		// 		if(min_abs_diff > component_diff_x_cur[i]) min_abs_diff = component_diff_x_cur[i];

		// 	if(min_abs_diff <= AdaptivePrecisionBound)
		// 		mat_increase_precision(A);

		// 	vector_copy(n, component_diff_x_cur, component_diff_x_prev);
		// }
		// update x vector
		// floatm_copy(n, x, x_prev);
		
		iter++;
		(*in_iter)++;
	}

    FREE(r);
    FREE(p);
    FREE(z);
    FREE(tr);
	// FREE(w);
	// FREE(x_prev); FREE(component_diff_x_prev);
}