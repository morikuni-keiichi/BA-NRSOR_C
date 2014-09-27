double drand(double l, double h);
double nrm2(double x[], int k);
void drotg(double *h1, double *h2, double *c, double *s);
void read_mat(int argc, char **argv);
void read_prm();
void output(int iter, double relres[], double t_tot, double x[]);
void BAGMRES(int *iter, double relres[], double x[]);
void NRSOR(double r[], double w[]);