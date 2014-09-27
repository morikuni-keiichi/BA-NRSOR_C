double drand(double l, double h);
double nrm2(double x[], int k);
void drotg(double *h1, double *h2, double *c, double *s);
void read_mat(int argc, char **argv);
void read_prm();
void output(int iter, double relres[], double t_tot, double x[]);
void BAGMRES(int *iter, double relres[], double x[]);
void NRSOR(double rhs[], double x[]);
void opNRSOR(double rhs[], double x[]);


extern double eps, omg, one, zero;
extern double *AC, *AR, *b, *Aei;
extern int *ia, *ip, *ja, *jp, m, maxit, n, nin, nnz;