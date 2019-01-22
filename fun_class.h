#ifndef fun_class
#define fun_class

int agaus(double a[],double b[], int n);
double *solve(double cor1[4],double cor2[4]);
double *A(double rad,double theta, double U);
double *solv(double rad, double theta);
double U(double rad,double theta);
double momen(double rad,double theta);
double *func(double rad,double theta,double psi,double u);
void NRSOLVE(double cor[4],double cor1[4]);
double momen_0(double rad,double theta);
double *gb_solve(double cor1[4],double cor2[4]);
double *dif_fun(double rad,double theta);
/*double *solv_2(double rad, double theta);*/
double *solv(double rad, double theta);

#endif