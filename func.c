#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern double one, zero;

double drand(double l, double h) {

	return ((double)rand() * (h-l)) / (double)RAND_MAX + l;
}


double nrm2(double *x, int k) {

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