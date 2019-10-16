#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <math.h>

double metro_x (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double *E, double dx, double beta, double *type);
double metro_p (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double *E, double dp, double beta, double *type);

#endif
