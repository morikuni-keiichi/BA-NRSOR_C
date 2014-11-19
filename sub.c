#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "header.h"

#define min(a, b) ((a) < (b) ? (a) : (b))


// Read parameter values
void read_prm() {

	FILE *file; 

	 // Check if prm.dat opens
  	if ((file = fopen("prm.dat", "r")) == NULL) {
  		exit(1);
  	}

  	fscanf(file, "%lf", &eps);
    fscanf(file, "%lf", &omg);
  	fscanf(file, "%d", &maxit);
    fscanf(file, "%d", &nin);
  	fclose(file);
    
}


void read_mat(int argc, char *argv[]) {

	int i, j, k;
	char *fi_AC = "AC.ccs", *fi_ia = "ia.ccs", *fi_jp = "jp.ccs";	
	char filename[256];
	FILE *file;

	if (argc > 2) {
		fprintf(stderr, "Error: Too many arguments.");
		exit(1);
	} else if (argc < 2) {
		fprintf(stderr, "Error: Too few arguments.");
		exit(1);
	}

	// Get AC.ccs
	strcpy(filename, argv[1]);
	strcat(filename, fi_AC);

	// Check if AC.ccs opens
	if ((file = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Failed to open AC.ccs");
		exit(1);
	} else {
		fscanf(file,"%d%d%d", &m, &n, &nnz);

		// Allocate AC
		if ((AC = (double *)malloc(sizeof(double) * (nnz))) == NULL) {
	    	fprintf(stderr, "Failed to allocate AC");
		  	exit(1);
		} else {
			for (k = 0; k<nnz; k++) fscanf(file, "%lf", &AC[k]);
			fclose(file);
		}
	}

	// Allocate ia
	if ((ia = (int *)malloc(sizeof(int) * (nnz))) == NULL) {
		fprintf(stderr, "Failed to allocate ia"); 
		exit(1);
	}
	
	// Allocate jp
	if ((jp = (int *)malloc(sizeof(int) * (n+1))) == NULL) {
		fprintf(stderr, "Failed to allocate jp");
		exit(1);
	}

	// Allocate Aei
	if ((Aei = (double *)malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "Failed to allocate Aei");
		exit(1);
	}

	if (maxit > min(m, n) || maxit < 1) {
		maxit = min(m, n);
	}

	// Get ia.crs
	strcpy(filename, argv[1]);
  	strcat(filename, fi_ia);

  	// Check if ia.ccs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stderr, "Failed to open ia.ccs");
  		exit(1);
  	} else {	// Load ia.ccs
  		for (k=0; k<nnz; k++) {
  			fscanf(file, "%d", &ia[k]);
  			ia[k] -= 1;
  		}
  		fclose(file); // Close ia.ccs
  	}  	

  	// Get filename jp.ccs
	strcpy(filename, argv[1]);
  	strcat(filename, fi_jp);

	// Check if jp.ccs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stderr, "Failed to open jp.ccs");
  		exit(1);
  	} else {	// Load jp.ccs
  		for (j=0; j<n+1; j++) {
  			fscanf(file, "%d", &jp[j]);
	  		jp[j] -= 1;
  		}  		
  		fclose(file); // Close jp.ccs
  	}	

	// Allocate b
	if ((b = (double *)malloc(sizeof(double) * (m))) == NULL) {
		fprintf(stderr, "Failed to allocate b");
		exit(1);
	} else { 
		for (i=0; i<m; i++) b[i] = drand(-one, one);
  	}

}


void output(int iter, double *relres, double t_tot, double *x) {

	double nrmATb, nrmATr, tmp, r[m], ATr[n];
	int i, j, k, k1, k2, l;
	FILE *file;	

	// A^T b
  	for (j=0; j<n; j++) {
		tmp = zero;
		k1 = jp[j];
		k2 = jp[j+1];
		for (l=k1; l<k2; l++) tmp += AC[l]*b[ia[l]];
		ATr[j] = tmp;
	}

  	nrmATb = nrm2(ATr, n);

	// A x
	for (i=0; i<m; i++) r[i] = zero; // Initialize
	for (j=0; j<n; j++) {
		tmp = x[j];
		for (l=jp[j]; l<jp[j+1]; l++) r[ia[l]] += tmp*AC[l];
 	}

  	// r = b - A x
	for (i=0; i<m; i++) r[i] = b[i] - r[i];

  	// ATr
	for (j=0; j<n; j++) {
		tmp = zero;
		k1 = jp[j];
		k2 = jp[j+1];
		for (l=k1; l<k2; l++) tmp += AC[l]*r[ia[l]];
		ATr[j] = tmp;
	}

  	nrmATr = nrm2(ATr, n);

  	printf("omega: %6.3f\n", omg);
	printf("# of outer iterations: %d\n", iter);
	printf("# of inner iterations: %d\n", nin);
	printf("CPU time: %16.5f\n", t_tot);
	printf("Relative residual: %18.16e\n", relres[iter-1]);
	printf("Actual relative residual: (ATr): %18.16e\n", nrmATr/nrmATb);

	// Check if info.dat opens
  	if ((file = fopen("info.dat", "w")) == NULL) {
  		fprintf(stderr, "Failed to open info.dat");
  		exit(1);
  	}

	fprintf(file, "omega: %6.3f\n", omg);
	fprintf(file, "# of outer iterations: %d\n", iter);
	fprintf(file, "# of inner iterations: %d\n", nin);
	fprintf(file, "CPU time: %16.5f\n", t_tot);
	fprintf(file, "Relative residual: %18.16e\n", relres[iter-1]); 
	fprintf(file, "Actual relative residual: (ATr): %18.16e\n", nrmATr/nrmATb);

  	fclose(file);

  	// Check if relres.dat opens
  	if ((file = fopen("reshis.dat", "w")) == NULL) {
  		fprintf(stderr, "Failed to open relres.dat");
  		exit(1);
  	} else { 
  		for (k=0; k<iter; k++) fprintf(file, "%.15e\n", relres[k]);
  		fclose(file);
  	}  	

  	// Check if solution.dat opens
  	if ((file = fopen("solution.dat", "w")) == NULL) {
  		fprintf(stderr, "Failed to open solution.dat");
  		exit(1);
  	} else {
  		for (j=0; j<n; j++) fprintf(file, "%.15e\n", x[j]);
  		fclose(file);
  	}

}