#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern double one, zero;

double drand(double l, double h) {

// Purpose
// =======
// 
// drand returns a random double variable
// 
// 
// 
// 
// =====================================================================

	return ((double)rand() * (h-l)) / (double)RAND_MAX + l;
}


double nrm2(double *x, int k) {

// Purpose
// =======
// 
// nrm2 returns the euclidean norm of a vector via the function
// name, so that
// 
// nrm2 := sqrt( x'*x )
// 
// Further Details
// ===============
// 
// -- This version written on 25-October-1982.
// Modified on 14-October-1993 to inline the call to DLASSQ.
// Sven Hammarling, Nag Ltd.
// translated into C and changed by Keiichi Morikun
// 
// References:
// 
// Jack Dongarra, Jim Bunch, Cleve Moler, and Pete Stewart,
// LINPACK User's Guide,
// Society for Industrial and Applied Mathematics (SIAM),  
// Philadelphia, 1979,
// ISBN-10: 089871172X,
// ISBN-13: 978-0-898711-72-1.
// 
// Charles Lawson, Richard Hanson, David Kincaid, and Fred Krogh,
// Algorithm 539, 
// Basic Linear Algebra Subprograms for Fortran Usage,
// ACM Transactions on Mathematical Software, 
// Volume 5, Number 3, September 1979, Pages 308-323.
// 
// =====================================================================

	double absxi, scale = zero, ssq = one, tmp;
	int i;

	for (i=0; i<k; i++) {
		if (x[i] != zero) {
	  		absxi = fabs(x[i]);
	    	if (scale <= absxi) {
	    		tmp = scale/absxi;
	    		ssq = one + ssq*tmp*tmp;
	    		scale = absxi;
	    	} else {
	    		tmp = absxi/scale;
	    		ssq = ssq + tmp*tmp;
			}
		}
	}

	return scale*sqrt(ssq);
}


void drotg(double *da, double *db, double *c, double *s)
{

// Purpose
// =======
// 
// drotg construct givens plane rotation.
// 
// Further Details
// ===============
// 
// jack dongarra, linpack, 3/11/78.
// translated into C and changed by Keiichi Morikuni
// 
// References:
// 
// Jack Dongarra, Jim Bunch, Cleve Moler, and Pete Stewart,
// LINPACK User's Guide,
// Society for Industrial and Applied Mathematics (SIAM),  
// Philadelphia, 1979,
// ISBN-10: 089871172X,
// ISBN-13: 978-0-898711-72-1.
// 
// Charles Lawson, Richard Hanson, David Kincaid, and Fred Krogh,
// Algorithm 539, 
// Basic Linear Algebra Subprograms for Fortran Usage,
// ACM Transactions on Mathematical Software, 
// Volume 5, Number 3, September 1979, Pages 308-323.
// 
// =====================================================================    

	double r, roe, scale, z;

	roe = *db;

    if (fabs(*da) > fabs(*db)) roe = *da;

    scale = fabs(*da) + fabs(*db);

    if (scale != zero) {
   	
	   	r = scale*sqrt(pow(*da/scale, 2.0) + pow(*db/scale, 2.0));
		 
		if (roe<0) r = -r;
	    *c = *da / r;
	    *s = *db / r;
	    z = one;

	    if (fabs(*da) > fabs(*db)) z = *s;

	    if (fabs(*db) >= fabs(*da) && *c != zero) z = one / *c;
		
		*da = r;
		*db = z;

 	} else {

 		*c = one;
    	*s = zero; 
    	r = zero;
    	z = zero;

    	*da = r;
    	*db = z;
    }

}