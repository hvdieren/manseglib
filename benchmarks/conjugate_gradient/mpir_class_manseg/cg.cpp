#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter)
{
    int iter = 0;
    FLOAT2 alpha, beta, rho, tau, tol;

	int initial_full_iter_limit = -1;

    FLOAT *r = ALLOC(FLOAT, n);     FLOAT *r_old = CALLOC(FLOAT, n);
    FLOAT *p = ALLOC(FLOAT, n);
    FLOAT *z = ALLOC(FLOAT, n);

    FLOAT *tr = ALLOC(FLOAT, n);

    floatm_mult(A, x, r); // r = Ax
    floatm_xpby(n, b, -1.0, r); // r0 = b - Ax
    floatm_copy(n, r, r_old);

    if (M) {
        floatm_mult(M, r, p);
        rho = floatm_dot(n, r, p);
        tol = floatm_norm2(n, r);
    } else { // small optimization
        floatm_copy(n, r, p);
        rho = floatm_dot(n, r, r);
        tol = sqrt(rho);
    }

    int step = 0;
    FLOAT2 residual = tol, prev_residual = tol;

    while ((iter < maxiter) && (tol > umbral)) {
        // alpha = (r,z) / (Ap,p)
        // for (int i = 0; i < n; i++) p[i] = double(p[i]);
        floatm_mult(A, p, z); // z = Ap

        // compute true residual
        if (step < step_check) step++;
        else {
            floatm_mult(A, x, tr); // tr = Ax
            floatm_xpby(n, b, -1.0, tr); // tr = b - Ax
            residual = floatm_norm2(n, tr);

            // FLOAT2 residual_diff = prev_residual - residual;
            // printf("res_def=%e\n", residual_diff);
            // // printf("sqrt(res_diff)=%e\n", sqrt(residual_diff));
            // prev_residual = residual;

            // if(!A->useTail && residual_diff < AdaptivePrecisionBound) {
            //     printf("---- switch pres at iter=%d , res=%e, res_diff=%e", *in_iter, residual, residual_diff);
            //     mat_switch(A);
            // }

            printf("# rescheck: total_cg_iter=%d current_iter=%d tol=%e resid=%e\n", *in_iter, iter, (double)tol, (double)residual);
            if (residual / tol > 10) {
                break;
            }
            step = 1;
        }

        alpha = rho / floatm_dot(n, z, p);
        
        // x = x + alpha * p
        floatm_axpy(n, alpha, p, x);
        
        // r = r - alpha * Ap
        floatm_axpy(n, -alpha, z, r);

        // this r may be what you're looking for
        // above is the last place where r is modified
        // so here you could check some form of rold/rnew
        // though, that's basically what the step-check residual is..
        // if(step > (step_check - 1)) {
        //     FLOAT2 res_diff = 1.0/sqrt(floatm_dot(n, r, r_old));
        //     printf("iter=%d, res_diff=%e\n", *in_iter, res_diff);
        //     floatm_copy(n, r, r_old);

        //     if(res_diff < AdaptivePrecisionBound) mat_switch(A);
        // }

        // apply preconditioner
        if (M) {
            floatm_mult(M, r, z);
            tau = floatm_dot(n, r, z);
            tol = floatm_norm2(n, r);
        } else {
            tau = floatm_dot(n, r, r);
            tol = sqrt(tau);
        }
        // beta = (r,z) / rho
        beta =  tau / rho;
        rho = tau;
        // p = z + beta * p
        if (M) floatm_xpby(n, z, beta, p);
        else floatm_xpby(n, r, beta, p);
        iter++;
        (*in_iter)++;
        /* printf ("%d %d %d %e %e\n", *in_iter, iter, bits, (double)tol, (double)residual);
        printf("alpha %e %e\n", alpha.value(), alpha.error());
        printf("beta %e %e\n", beta.value(), beta.error());
        printf("rho %e %e\n", rho.value(), rho.error());
        printf("tau %e %e\n", tau.value(), tau.error()); */

		if(*in_iter == initial_full_iter_limit) {
			printf("reduced precision after %d iterations\n", *in_iter);
			mat_reduce_precision(A);
		}
    }

    FREE(r);            FREE(r_old);
    FREE(p);
    FREE(z);
    FREE(tr);
}