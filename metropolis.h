#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <math.h>

double aceptacion (double *v, double *v_new, double *E, double dE, double beta, int N, int target);
double metro_x (double *x, double *p, double *E, double *tabla_V_LJ, double *tabla_V_P, double rc2, double sc2, double dr2, double ds2, double q0, double p0, double L, double beta, double dx, int N);
double metro_p (double *x, double *p, double *E, double *tabla_V_LJ, double *tabla_V_P, double sc2, double ds2, double q0, double p0, double L, double beta, double dp, int N);

#endif
