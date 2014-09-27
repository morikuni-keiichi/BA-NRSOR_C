#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "header.h"

extern double eps, one, omg, zero;
extern double *AC, *AR, *b, *Aei;
extern int *ia, *ip, *ja, *jp, m, maxit, n, nin, nnz;

void BAGMRES(int *iter, double *relres, double *x) {

	double **H = NULL, **V = NULL;
	double c[maxit], g[maxit+1], r[m], s[maxit], w[n], y[maxit];
	double beta, nrmATb, nrmATr, tmp, t_tot, Tol;
	int i, j, k, k1, k2, l;

	// Allocate H[maxit]
	if ((H = malloc(sizeof(double) * (maxit))) == NULL) {
		fprintf(stderr, "Failed to allocate H");
		exit(1);
	}

	// Allocate H[maxit][maxit+1]
	for (k=0; k<maxit; k++) {
		if ((H[k] = (double *)malloc(sizeof(double) * (k+2))) == NULL) {
			fprintf(stdout, "Failed to allocate H");
			exit(1);
		}
	}

	// Allocate V[maxit+1]
	if ((V = malloc(sizeof(double) * (maxit+1))) == NULL) {
		fprintf(stdout, "Failed to allocate V");
		exit(1);
	}

	// Allocate V[maxit+1][n]
	for (j=0; j<maxit+1; j++) {
		if ((V[j] = (double *)malloc(sizeof(double) * n)) == NULL) {
			fprintf(stdout, "Failed to allocate V");
			exit(1);;
		}
	}
	
	// Norm of the jth column of A for NR-SOR
	for (j=0; j<n; j++) Aei[j] = zero;
	for (i=0; i<m; i++) {
		for (l=ip[i]; l<ip[i+1]; l++) {
			j = ja[l];
			Aei[j] = Aei[j] + AR[l]*AR[l];
		}
	}

	// Take the inverse of ||aj||
	for (j=0; j<n; j++) {
		if (Aei[j] > zero) {
			Aei[j] = omg / Aei[j];
		} else {
			fprintf(stderr, "warning: ||aj|| = 0.0");
			exit(1);
		}
	}

	for (j=0; j<n; j++) w[j] = zero; // Initilize
	// w = A^T b
	for (i=0; i<m; i++) {
  		tmp = b[i];
  		for (l=ip[i]; l<ip[i+1]; l++) {
  			j = ja[l];
  			w[j] = w[j] + tmp*AR[l];
  		}
  	}

  	// norm of A^T b
  	nrmATb = nrm2(w, n);

  	// Stopping criterion
  	Tol = eps * nrmATb;

  	// r = b  (x0 = 0)
  	for (i=0; i<m; i++) r[i] = b[i];

  	// NR-SOR inner iterations: w = B r
  	NRSOR(r, w);

  	// beta = norm(Bb)
  	beta = nrm2(w, n);

  	// Normalize
  	tmp = one / beta;
  	for (j=0; j<n; j++) V[0][j] = tmp * w[j];

  	// beta e1
  	g[0] = beta;

  	// Main loop
  	for (k=0; k<maxit; k++) {

		// r = A v 	
  		for (i=0; i<m; i++) {
			tmp = zero;
			k1 = ip[i];
			k2 = ip[i+1];
			for (l=k1; l<k2; l++) tmp += AR[l]*V[k][ja[l]];
			r[i] = tmp;
		}

		// NR-SOR inner iterations: w = B r
		NRSOR(r, w);

		// Modified Gram-Schmidt orthogonzlization
		for (i=0; i<k+1; i++) {
			tmp = zero;
			for (j=0; j<n; j++) tmp += w[j]*V[i][j];
			H[k][i] = tmp;
			for (j=0; j<n; j++) w[j] = w[j] - tmp*V[i][j];
		}

		// h_{kL1, k}
		H[k][k+1] = nrm2(w, n);

		// Check breakdown
		if (H[k][k+1] > zero) {
			tmp = one / H[k][k+1];
			for (j=0; j<n; j++) V[k+1][j] = tmp * w[j];
		} else {
			fprintf(stderr, "Breakdown at%dth step\n", k+1);
			fprintf(stderr, "h_k+1, k = %.15e\n", H[k][k+1]);
			fflush(stderr);
		}

		// Apply Givens rotations
		for (i=0; i<k; i++) {
			tmp = c[i]*H[k][i] + s[i]*H[k][i+1];
			H[k][i+1] = -s[i]*H[k][i] + c[i]*H[k][i+1];
			H[k][i] = tmp;
		}	

		// Compute Givens rotations
		drotg(&H[k][k], &H[k][k+1], &c[k], &s[k]);

		// Apply Givens rotations
		g[k+1] = -s[k] * g[k];
		g[k] = c[k] * g[k];

		nrmATr = fabs(g[k+1]);
		relres[k] = nrmATr;

		if (nrmATr < Tol || k == maxit-1) {

			y[k] = g[k] / H[k][k];
			for (i=k-1; i>-1; i--) {				
				tmp = zero;
				for (l=i+1; l<k; l++) tmp += H[l][i]*y[l];				
				y[i] = (g[i] - tmp) / H[i][i];
			}

			for (j=0; j<n; j++) x[j] = zero;
			for (l=0; l<k+1; l++) {
				for (j=0; j<n; j++) x[j] += V[l][j]*y[l];
			}

			// r = A x 	
	  		for (i=0; i<m; i++) {
				tmp = zero;
				k1 = ip[i];
				k2 = ip[i+1];
				for (l=k1; l<k2; l++) tmp += AR[l]*x[ja[l]];
				r[i] = tmp;
			}
			
			// w = A^T r
			for (i=0; i<m; i++) r[i] = b[i] - r[i];

			for (j=0; j<n; j++) w[j] = zero; // Initilize
			for (i=0; i<m; i++) {
				tmp = r[i];
				for (l=ip[i]; l<ip[i+1]; l++) {
					j = ja[l];
					w[j] = w[j] + tmp*AR[l];
				}
		 	}		 	

		  	// printf("%.15e\n", nrm2(w, n)/nrmATb);

		  	if (nrm2(w, n) < Tol || k == maxit-1) {
				*iter = k+1;

				for (j=0; j<maxit+1; j++) {
					free(V[j]);
				}		
				free(V);

			  	for (k=0; k<maxit; k++) {
					free(H[k]);
				}
				free(H);

				return;
			}
		}

	}

}


void NRSOR(double *rhs, double *x)
{
	double d;
	int i, j, k, k1, k2, l;

	for (j=0; j<n; j++) x[j] = zero;

	for (k=0; k<nin; k++) {
		for (j=0; j<n; j++) {
			k1 = jp[j];
			k2 = jp[j+1];
			d = zero;
			for (l=k1; l<k2; l++) d += AC[l]*rhs[ia[l]];
			d = d * Aei[j];			
			x[j] = x[j] + d;
			for (l=k1; l<k2; l++) {
				i = ia[l];
				rhs[i] = rhs[i] - d*AC[l];
			}
		}
	}
}