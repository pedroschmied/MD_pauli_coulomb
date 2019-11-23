#ifndef INTERACCION_H
#define INTERACCION_H
#include <math.h>

int interaccion_PCB(double *delta_r, double L);
int pos_PBC(double *x, int i, double L);
double hamiltoneano (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double m);
double eval_LJ (double *tabla_V_LJ, double rij2, double rc2, double dr2);
double eval_P (double *tabla_V_P, double rij2, double pij2, double q0, double p0, double sc2, double ds2);
double  variacion_E_x (double *x, double *x_new, double *p, double *tabla_V_LJ, double *tabla_V_P, double rc2, double sc2, double dr2, double ds2, double q0, double p0, double L, int N, int target);
double  variacion_E_p (double *x, double *p, double *p_new, double *tabla_V_P, double sc2, double ds2, double q0, double p0, double L, int N, int target, double m);

#endif
