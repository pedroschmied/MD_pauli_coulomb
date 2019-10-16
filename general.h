#ifndef GENERAL_H
#define GENERAL_H
#include <math.h>

double aleatorio();
double mean (double *v, int o, int n, int k);
double mean2 (double *v, int o, int n, int k);
double std2 (double *v, int o, int n, int k);
int delta_x(double *x, int i, double *y, int j, double *delta_r);
double gaussiana(double mu, double sigma);
double norma2(double *x);


#endif
