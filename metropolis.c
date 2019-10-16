#include "interaccion.h"
#include "metropolis.h"

double metro_x (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double *E, double dx, double beta, double *type)
{
	double *x_new;
	x_new = (double*) malloc(3 * N * sizeof(double));
	int z, target;
	double a = aleatorio();
	target = (int)(a * ((double)N - 1.0)); //eligo la partícula
	for (z = 0; z < 3 * N; z++)
	{
		*(x_new + z) = *(x + z);
	}
	for (z = 0; z < 3; z++)
	{
		a = aleatorio();
		*(x_new + 3 * target + z) = (1.0 - 2.0 * a) * dx + *(x + 3 * target + z);
	}
	pos_PBC(x_new, target, N, L);
	double dE;
	dE = variacion_E  (x, p, x_new, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, target, type);
	double aceptacion = 0.0;
	if (dE <= 0.0)
	{
		for (z = 0; z < 3; z++)
		{
			*(x + 3 * target + z) = *(x_new + 3 * target + z);
		}
		*E += dE;
		aceptacion = 1.0;
	}
	else
	{
		a = aleatorio();
		if(a < exp(-dE * (double) N * beta))
		{
			for (z = 0; z < 3; z++)
			{
				*(x + 3 * target + z) = *(x_new + 3 * target + z);
			}
			*E += dE;
			aceptacion = 1.0;
		}
	}
	free(x_new);
	return aceptacion;
}

double metro_p (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double *E, double dp, double beta, double *type)
{
	double *p_new;
	p_new = (double*) malloc(3 * N * sizeof(double));
	int z, target;
	double a = aleatorio();
	target = (int)(a * ((double)N - 1.0)); //eligo la partícula
	for (z = 0; z < 3 * N; z++)
	{
		*(p_new + z) = *(p + z);
	}
	for (z = 0; z < 3; z++)
	{
		a = aleatorio();
		*(p_new + 3 * target + z) = (1.0 - 2.0 * a) * dp + *(p + 3 * target + z);
	}
	double dE;
	dE = variacion_E  (x, p, x, p_new, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, target, type);

	double aceptacion = 0.0;
	if (dE <= 0.0)
	{
		for (z = 0; z < 3; z++)
		{
			*(p + 3 * target + z) = *(p_new + 3 * target + z);
		}
		*E += dE;
		aceptacion = 1.0;
	}
	else
	{
		a = aleatorio();
		if(a < exp(-dE * (double) N * beta))
		{
			for (z = 0; z < 3; z++)
			{
				*(p + 3 * target + z) = *(p_new + 3 * target + z);
			}
			*E += dE;
			aceptacion = 1.0;
		}
	}
	free(p_new);
	return aceptacion;
}
