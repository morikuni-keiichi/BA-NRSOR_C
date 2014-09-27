#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "header.h"

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
			fprintf(stderr, "Failed to allocate H");
			exit(1);
		}
	}

	// Allocate V[maxit+1]
	if ((V = malloc(sizeof(double) * (maxit+1))) == NULL) {
		fprintf(stderr, "Failed to allocate V");
		exit(1);
	}

	// Allocate V[maxit+1][n]
	for (j=0; j<maxit+1; j++) {
		if ((V[j] = (double *)malloc(sizeof(double) * n)) == NULL) {
			fprintf(stderr, "Failed to allocate V");
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
			Aei[j] = one / Aei[j];
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
  	// NRSOR(r, w);
  	opNRSOR(r, w); 
  
  	for (j=0; j<n; j++) Aei[j] = omg * Aei[j];

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

		relres[k] = fabs(g[k+1]) / beta; 

		if (relres[k] < eps || k == maxit-1) {

			// Derivation of the approximate solution x_k
			// Backward substitution		
			y[k] = g[k] / H[k][k];
			for (i=k-1; i>-1; i--) {				
				tmp = zero;
				for (l=i+1; l<k; l++) tmp += H[l][i]*y[l];				
				y[i] = (g[i] - tmp) / H[i][i];
			}

			// x = V y
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

		 	nrmATr = nrm2(w, n);

		 	// Convergence check
		  	if (nrmATr < Tol || k == maxit-1) {
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

	for (k=1; k<=nin; k++) {
		for (j=0; j<n; j++) {
			k1 = jp[j];
			k2 = jp[j+1];
			d = zero;
			for (l=k1; l<k2; l++) d += AC[l]*rhs[ia[l]];
			d = d * Aei[j];			
			x[j] = x[j] + d;
			if (k == nin && j == n-1) return;
			for (l=k1; l<k2; l++) {
				i = ia[l];
				rhs[i] = rhs[i] - d*AC[l];
			}
		}
	}
}


void opNRSOR(double *rhs, double *x)
{
	double d, e, res1, res2 = zero, tmp, y[n], tmprhs[m];
	int i, ii, j, k, k1, k2, l;

	// Initilize
	for (i=0; i<m; i++) tmprhs[i] = rhs[i];

	for (j=0; j<n; j++) {
		x[j] = zero;
		y[j] = zero;	
	}

	// Tune the number of inner iterations 
	for (k=1; k<=100; k++) {

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

		d = zero;
		for (j=0; j<n; j++) { 
			tmp = fabs(x[j]);
			if (d < tmp) d = tmp;
		}

		e = zero;
		for (j=0; j<n; j++) { 
			tmp = fabs(x[j] - y[j]);
			if (e < tmp) e = tmp;
		}

		if (e<1.0e-1*d || k == 100) {
			nin = k;
			break;

		}

		for (j=0; j<n; j++) y[j] = x[j];

	}

	// Tune the relaxation parameter
	for (k=19; k>0; k--) {
		omg = 1.0e-1 * (double)(k); // omg = 1.9, 1.8, ..., 0.1

		for (i=0; i<m; i++) rhs[i] = tmprhs[i];

		for (j=0; j<n; j++) x[j] = zero;

		for (i=1; i<=nin; i++) {
			for (j=0; j<n; j++) {
				k1 = jp[j];
				k2 = jp[j+1];
				d = zero;
				for (l=k1; l<k2; l++) d += AC[l]*rhs[ia[l]];
				d = omg * d * Aei[j];			
				x[j] = x[j] + d;
				for (l=k1; l<k2; l++) {
					ii = ia[l];
					rhs[ii] = rhs[ii] - d*AC[l];
				}
			}
		}

		res1 = nrm2(rhs, m);

		if (k < 19) {
			if (res1 > res2) {
				omg = omg + 1.0e-1;
				for (j=0; j<n; j++) x[j] = y[j];					
				return;
			} else if (k == 1) {
				omg = 1.0e-1;
				return;
			}
		}

		res2 = res1;

		for (j=0; j<n; j++) y[j] = x[j];
	}
}