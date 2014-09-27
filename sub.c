#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "header.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

extern double eps, omg, one, zero;
extern double *AC, *AR, *b, *Aei;
extern int *ia, *ip, *ja, *jp, m, maxit, n, nin, nnz;

void read_prm() {

	FILE *file; 

	 // Chech if AR.crs opens
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

	int i, j, k, m_AC, m_AR, n_AC, n_AR, nnz_AC, nnz_AR;
	char filename[256];
	char *fi_AC = "AC.ccs", *fi_ia = "ia.ccs", *fi_jp = "jp.ccs";
	char *fi_AR = "AR.crs", *fi_ip = "ip.crs", *fi_ja = "ja.crs";
	FILE *file;

	if (argc > 2) {
		fprintf(stderr, "Error: Too many arguments.");
		exit(1);
	} else if (argc < 2) {
		fprintf(stderr, "Error: Too few arguments.");
		exit(1);
	}

	strcpy(filename, argv[1]);
  	strcat(filename, fi_AC);

  	// Check if AC.ccs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stdout, "Failed to open AC.ccs");
  		exit(1);
  	}

  	fscanf(file,"%d%d%d", &m_AC, &n_AC, &nnz_AC);

  	// Allocation
  	if ((AC = (double *)malloc(sizeof(double) * (nnz_AC))) == NULL) {
		fprintf(stdout, "Failed to allocate AC");
		exit(1);
	}

	for (k = 0; k<nnz_AC; k++) fscanf(file, "%lf", &AC[k]);

  	fclose(file);

	strcpy(filename, argv[1]);
  	strcat(filename, fi_AR);

  	// Chech if AR.crs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stdout, "Failed to open AR.crs");
  		exit(1);
  	}

  	fscanf(file,"%d%d%d", &m_AR, &n_AR, &nnz_AR);

  	if (m_AR != m_AC || n_AR != n_AC || nnz_AR != nnz_AC) {
  		fprintf(stderr, "CRS and CCS files are not consistent.");
  		exit(1);
	} else {
		m = m_AR;
		n = n_AR;
		nnz = nnz_AR;
	}

  	// Allocation
  	if ((AR = (double *)malloc(sizeof(double) * (nnz))) == NULL) {
		fprintf(stdout, "Failed to allocate AR");
		exit(1);
	}

	// Allocation
	if ((ia = (int *)malloc(sizeof(int) * (nnz))) == NULL) {
		fprintf(stderr, "Failed to allocate ia"); 
		exit(1);
	}

	// Allocation
	if ((ip = (int *)malloc(sizeof(int) * (m+1))) == NULL) {
		fprintf(stderr, "Failed to allocate ip"); 
		exit(1);
	}

	// Allocation
	if ((ja = (int *)malloc(sizeof(int) * (nnz))) == NULL) {
		fprintf(stderr, "Failed to allocate ja");
		exit(1);
	}

	// Allocation
	if ((jp = (int *)malloc(sizeof(int) * (n+1))) == NULL) {
		fprintf(stderr, "Failed to allocate jp");
		exit(1);
	}

	// Allocation
	if ((Aei = (double *)malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "Failed to allocate Aei");
		exit(1);
	}

	if (maxit > min(m, n) || maxit < 1) {
		maxit = min(m, n);
	}

	for (k = 0; k<nnz; k++) fscanf(file, "%lf", &AR[k]);

	fclose(file);

	// Get ia.crs
	strcpy(filename, argv[1]);
  	strcat(filename, fi_ia);

  	// Chech if ia.ccs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stdout, "Failed to open ia.ccs");
  		exit(1);
  	}

  	// Load ia.ccs
  	for (k=0; k<nnz; k++) {
  		fscanf(file, "%d", &ia[k]);
  		ia[k] = ia[k] - 1;
  	}

  	// Close ia.ccs
  	fclose(file);

	// Get ip.crs
	strcpy(filename, argv[1]);
  	strcat(filename, fi_ip);

  	// Chech if ip.crs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stdout, "Failed to open ip.crs");
  		exit(1);
  	}

  	// Load ip.crs
  	for (i = 0; i<m+1; i++) {
  		fscanf(file, "%d", &ip[i]);
  		ip[i] = ip[i] - 1;
  	}

  	// Close ip.crs
  	fclose(file);

  	// Get filename ja.crs
	strcpy(filename, argv[1]);
  	strcat(filename, fi_ja);

	// Check if ja.crs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stdout, "Failed to open ja.crs");
  		exit(1);
  	}

  	// Load ja.crs
  	for (k = 0; k<nnz; k++) {
  		fscanf(file, "%d", &ja[k]);
  		ja[k] = ja[k] - 1;
  	}

	// Close ja.ccs
  	fclose(file);

  	// Get filename jp.ccs
	strcpy(filename, argv[1]);
  	strcat(filename, fi_jp);

	// Check if jp.ccs opens
  	if ((file = fopen(filename, "r")) == NULL) {
  		fprintf(stdout, "Failed to open jp.ccs");
  		exit(1);
  	}

  	// Load jp.ccs
  	for (j=0; j<n+1; j++) {
  		fscanf(file, "%d", &jp[j]);
  		jp[j] = jp[j] - 1;
  	}

	// Close jp.ccs
  	fclose(file);

	// Allocation
	if ((b = (double *)malloc(sizeof(double) * (m))) == NULL) {
		fprintf(stdout, "Failed to allocate b");
		exit(1);
	}

  	for (i=0; i<m; i++) {
  		b[i] = drand(-one, one);
  	}

}


void output(int iter, double *relres, double t_tot, double *x) {

	double nrmATb, nrmATr, tmp, r[m], ATr[n];
	int i, j, k, k1, k2, l;
	FILE *file;

	for (j=0; j<n; j++) ATr[j] = zero;
	for (i=0; i<m; i++) {
		tmp = b[i];
		for (l=ip[i]; l<ip[i+1]; l++) {
			j = ja[l];
			ATr[j] = ATr[j] + tmp*AR[l];
		}
  }		

  nrmATb = nrm2(ATr, n);

	// r = A*x
	for (i=0; i<m; i++) {
		tmp = zero;
		k1 = ip[i];
		k2 = ip[i+1];
		for (l=k1; l<k2; l++) tmp += AR[l]*x[ja[l]];
		r[i] = tmp;
	}

	for (i=0; i<m; i++) r[i] = b[i] - r[i];

	for (j=0; j<n; j++) ATr[j] = zero;
	for (i=0; i<m; i++) {
		tmp = r[i];
		for (l=ip[i]; l<ip[i+1]; l++) {
			j = ja[l];
			ATr[j] = ATr[j] + tmp*AR[l];
		}
  }	

  	nrmATr = nrm2(ATr, n);

  	tmp = one / nrmATb;
  	for (k=0; k<iter; k++) relres[k] = tmp * relres[k];

	printf("# of iterations: %d\n", iter);
	printf("CPU time: %7.5f\n", t_tot);
	printf("Relative residual: %18.16e\n", relres[iter-1]);
	printf("Actual relative residual: (ATr): %18.16e\n", nrmATr/nrmATb);

	// Chech if info.dat opens
  	if ((file = fopen("info.dat", "w")) == NULL) {
  		fprintf(stdout, "Failed to open info.dat");
  		exit(1);
  	}

  	fprintf(file, "# of iterations: %d\n", iter);
  	fprintf(file, "CPU time: %7.5f\n", t_tot);
	fprintf(file, "Relative residual: %18.16e\n", relres[iter-1]); 
	fprintf(file, "Actual relative residual: (ATr): %18.16e\n", nrmATr/nrmATb);

  	fclose(file);

  	// Chech if relres.dat opens
  	if ((file = fopen("reshis.dat", "w")) == NULL) {
  		fprintf(stderr, "Failed to open relres.dat");
  		exit(1);
  	}

  	for (k=0; k<iter; k++) fprintf(file, "%.15e\n", relres[k]);

  	fclose(file);

  	// Chech if solution.dat opens
  	if ((file = fopen("solution.dat", "w")) == NULL) {
  		fprintf(stderr, "Failed to open solution.dat");
  		exit(1);
  	}

  	for (j=0; j<n; j++) fprintf(file, "%.15e\n", x[j]);

  	fclose(file);

}