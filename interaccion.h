#ifndef INTERACCION_H
#define INTERACCION_H
#include <math.h>

int interaccion_PCB(double *delta_r, double L);
int pos_PBC(double *x, int i, int N, double L);
double hamiltoneano (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double *type);
double eval_interaccion (double *V_LJ, double *V_P, double *tabla_V_LJ, double *tabla_V_P, double rij2, double pij2, double p0, double q0, double rc2, double sc2, double ds2, double dr2, double *type, int i, int j);
double  variacion_E  (double *x, double *p, double *x_new, double *p_new, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, int target, double *type);

#endif
