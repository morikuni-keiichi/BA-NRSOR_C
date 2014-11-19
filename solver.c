#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "header.h"


// BA-GMRES
void BAGMRES(int *iter, double *relres, double *x) {

	double **H, **V, *c, *g, *r, *s, *w, *y;	
	double beta, nrmATb, nrmATr, tmp, t_tot, Tol;
	int i, j, k, k1, k2, l;

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
	
	// Allocate r
	if ((r = (double *)malloc(sizeof(double) * (m))) == NULL) {
		fprintf(stderr, "Failed to allocate r");
	}

	// Allocate w
	if ((w = (double *)malloc(sizeof(double) * (n))) == NULL) {
		fprintf(stderr, "Failed to allocate w");
	}

	// Allocate Aei
	if ((Aei = (double *)malloc(sizeof(double) * (n))) == NULL) {
		fprintf(stderr, "Failed to allocate Aei");
	}

		// Allocate c
	if ((c = (double *)malloc(sizeof(double) * (maxit))) == NULL) {
		fprintf(stderr, "Failed to allocate c");
	}

	// Allocate g
	if ((g = (double *)malloc(sizeof(double) * (maxit+1))) == NULL) {
		fprintf(stderr, "Failed to allocate g");
	}

	// Allocate s
	if ((s = (double *)malloc(sizeof(double) * (maxit))) == NULL) {
		fprintf(stderr, "Failed to allocate s");
	}

	// Allocate y
	if ((y = (double *)malloc(sizeof(double) * (maxit))) == NULL) {
		fprintf(stderr, "Failed to allocate y");
	}

	// Take the inverse of the norm of the jth column of A for NR-SOR
	for (j=0; j<n; j++) {		
		tmp = zero;
		k1 = jp[j];
		k2 = jp[j+1];	
		for (l=k1; l<k2; l++) tmp += AC[l]*AC[l];
		if (tmp > zero) {
			Aei[j] = one / tmp;
		} else {
			fprintf(stderr, "warning: ||aj|| = 0.0");
			exit(1);
		}
	}
	
	// // w = A^T b	
	for (j=0; j<n; j++) {
		tmp = zero;
		k1 = jp[j];
		k2 = jp[j+1];
		for (l=k1; l<k2; l++) tmp += AC[l]*b[ia[l]];
		w[j] = tmp;
	}

  	// norm of A^T b
  	nrmATb = nrm2(w, n);

  	// Stopping criterion
  	Tol = eps * nrmATb;

  	// r = b  (x0 = 0)
  	for (i=0; i<m; i++) r[i] = b[i];

  	// NR-SOR inner iterations: w = B r
  	fprintf(stdout, "Automatic NR-SOR inner-iteration parameter tuning\n");
  	fflush(stdout);
  	opNRSOR(r, w); 
  	fprintf(stdout, "Tuned\n");
  	fflush(stdout);
  
  	for (j=0; j<n; j++) Aei[j] *= omg;

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
		for (i=0; i<m; i++) r[i] = zero; // Initialize
		for (j=0; j<n; j++) {
			tmp = V[k][j];
			k1 = jp[j];
			k2 = jp[j+1];
			for (l=k1; l<k2; l++) r[ia[l]] += tmp*AC[l];
	 	}

		// NR-SOR inner iterations: w = B r
		NRSOR(r, w);

		// Modified Gram-Schmidt orthogonalization
		for (i=0; i<k+1; i++) {
			tmp = zero;
			for (j=0; j<n; j++) tmp += w[j]*V[i][j];			
			for (j=0; j<n; j++) w[j] -= tmp*V[i][j];
			H[k][i] = tmp;
		}

		// h_{kL1, k}
		tmp = nrm2(w, n);
		H[k][k+1] = tmp;

		// Check breakdown
		if (H[k][k+1] > zero) {
			tmp = one / tmp;
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

		H[k][k] = one / H[k][k];

		// Apply Givens rotations
		g[k+1] = -s[k] * g[k];
		g[k] *= c[k];

		relres[k] = fabs(g[k+1]) / beta; 

		if (relres[k] < eps) {

			// Derivation of the approximate solution x_k
			// Backward substitution		
			y[k] = g[k] * H[k][k];
			i = k;
			while (i--) {				
				tmp = zero;
				for (l=i+1; l<k; l++) tmp += H[l][i]*y[l];				
				y[i] = (g[i] - tmp) * H[i][i];
			}

			// x = V y
			for (j=0; j<n; j++) x[j] = zero;
			for (l=0; l<k; l++) {
				for (j=0; j<n; j++) x[j] += V[l][j]*y[l];
			}

			// r = A x
			for (i=0; i<m; i++) r[i] = zero; // Initialize
			for (j=0; j<n; j++) {
				tmp = x[j];
				for (l=jp[j]; l<jp[j+1]; l++) r[ia[l]] += tmp*AC[l];				
		 	}
			
			// r = b - Ax
			for (i=0; i<m; i++) r[i] = b[i] - r[i];

			// w = A^T r
		 	for (j=0; j<n; j++) {
				tmp = zero;
				k1 = jp[j];
				k2 = jp[j+1];
				for (l=k1; l<k2; l++) tmp += AC[l]*r[ia[l]];
				w[j] = tmp;
			}

		 	nrmATr = nrm2(w, n);

		 	// Convergence check
		  	if (nrmATr < Tol) {
				*iter = k+1;

				free(y);
				free(s);
				free(c);
				free(g);
				free(Aei);
				free(w);
				free(r);		

				k = maxit;
			  	while (k--) free(H[k]);
				free(H);

				j = maxit + 1;
				while (j--) free(V[j]);
				free(V);	

  				// printf("Required number of iterations: %d\n", (int)(*iter));
  				printf("Successfully converged.\n");  	

				return;
			}
		}

	}

	printf("Failed to converge.\n");

	k = k - 1;

	// Derivation of the approximate solution x_k
	// Backward substitution	

	y[k] = g[k] * H[k][k];
	i = k;
	while (i--) {	
		tmp = zero;
		for (l=i+1; l<k; l++) tmp += H[l][i]*y[l];				
		y[i] = (g[i] - tmp) * H[i][i];
	}

	// x = V y
	for (j=0; j<n; j++) x[j] = zero;
	for (l=0; l<k; l++) {
		for (j=0; j<n; j++) x[j] += V[l][j]*y[l];
	}

	*iter = k+1;

	free(y);
	free(s);
	free(c);
	free(g);
	free(Aei);
	free(w);
	free(r);		

	k = maxit;
  	while (k--) free(H[k]);
	free(H);

	j = maxit + 1;
	while (j--) free(V[j]);
	free(V);		

	return;

}


// NR-SOR inner-iteration preconditioning
void NRSOR(double *rhs, double *x)
{
	double d;
	int j, k, k1, k2, l;

	for (j=0; j<n; j++) x[j] = zero;

	k = nin;
	while (k--) {	
		for (j=0; j<n; j++) {
			k1 = jp[j];
			k2 = jp[j+1];
			d = zero;
			for (l=k1; l<k2; l++) d += AC[l]*rhs[ia[l]];
			d *= Aei[j];			
			x[j] += d;
			if (k == 0 && j == n-1) return;
			for (l=k1; l<k2; l++) rhs[ia[l]] -= d*AC[l];			
		}
	}
}


// Automatic parameter tuning for the NR-SOR inner-iteration preconditioning
void opNRSOR(double *rhs, double *x)
{
	double d, e, res1, res2 = zero, tmp, y[n], tmprhs[m];
	int i, ii, j, k, k1, k2, l;

	// Initialize
	for (i=0; i<m; i++) tmprhs[i] = rhs[i];

	for (j=0; j<n; j++) {
		x[j] = zero;
		y[j] = zero;	
	}

	// Tune the number of inner iterations 
	for (k=1; k<=50; k++) {

		for (j=0; j<n; j++) {
			k1 = jp[j];
			k2 = jp[j+1];
			d = zero;
			for (l=k1; l<k2; l++) d += AC[l]*rhs[ia[l]];
			d *= Aei[j];			
			x[j] += d;
			for (l=k1; l<k2; l++) rhs[ia[l]] -= d*AC[l];			
		}

		d = zero;
		e = zero;
		for (j=0; j<n; j++) { 
			tmp = fabs(x[j]);
			if (d < tmp) d = tmp;		
			tmp = fabs(x[j] - y[j]);
			if (e < tmp) e = tmp;
		}

		if (e<1.0e-1*d || k == 50) {
			nin = k;
			break;

		}

		for (j=0; j<n; j++) y[j] = x[j];

	}

	// Tune the relaxation parameter
	k = 20;
	while (k--)  {
		omg = 1.0e-1 * (double)(k); // omg = 1.9, 1.8, ..., 0.1

		for (i=0; i<m; i++) rhs[i] = tmprhs[i];

		for (j=0; j<n; j++) x[j] = zero;

		for (i=1; i<=nin; i++) {
			for (j=0; j<n; j++) {
				k1 = jp[j];
				k2 = jp[j+1];
				d = zero;
				for (l=k1; l<k2; l++) d += AC[l]*rhs[ia[l]];
				d *= omg * Aei[j];			
				x[j] += d;
				for (l=k1; l<k2; l++) rhs[ia[l]] -= d*AC[l];				
			}
		}

		res1 = nrm2(rhs, m);

		if (k < 19) {
			if (res1 > res2) {
				omg += 1.0e-1;
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