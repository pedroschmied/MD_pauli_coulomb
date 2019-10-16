#ifndef TABLAS_H
#define TABLAS_H
#include <math.h>

double tablas (double *tabla_V_LJ, double *tabla_V_P, int largo_tabla, double rc2, double sc2, double q0, double p0);
double interpol (double *tabla, double rij, double dr);
#endif
