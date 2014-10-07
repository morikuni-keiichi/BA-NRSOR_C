#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "header.h"

// Global variables
double eps, omg, one = 1.0, zero = 0.0;
double *AC = NULL, *b = NULL, *Aei = NULL;
int *ia = NULL, *jp = NULL, m, maxit, n, nin, nnz;

int main(int argc, char *argv[])
{

	double *relres = NULL, *x = NULL;
	double t_tot;	
	int iter;
	clock_t t1, t2;
	
	// Read parameters
	read_prm();

	// Read matrix data
	read_mat(argc, argv);

	// Allocate x
	if ((x = (double *)malloc(sizeof(double) * (n))) == NULL) {
		fprintf(stderr, "Failed to allocate x");
		exit(1);
	}

	// Allocate relres
	if ((relres = (double *)malloc(sizeof(double) * (maxit))) == NULL) {
		fprintf(stderr, "Failed to allocate relres");
		exit(1);
	}	
	
	// BA-GMRES
	t1 = clock();
  	BAGMRES(&iter, relres, x);
  	t2 = clock();
  	t_tot = (double)(t2 - t1) / CLOCKS_PER_SEC;

  	// Output results
  	output(iter, relres, t_tot, x);

  	// Deallocate
	free(AC);
	free(jp);
	free(ia);
	free(x);
	free(b);
	free(relres);

	return 0;
}