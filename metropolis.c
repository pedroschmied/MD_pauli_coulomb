#include "interaccion.h"
#include "metropolis.h"

double aceptacion (double *v, double *v_new, double *E, double dE, double beta, int N, int target)
{
	double a, si = 0.0;
	int z;
	if (dE <= 0.0)
	{
		for (z = 0; z < 3; z++)
		{
			*(v + 3 * target + z) = *(v_new + 3 * target + z);
		}
		*E += dE;
		si = 1.0;
	}
	else
	{
		a = aleatorio();
		if(a < exp(-dE * beta))
		{
			for (z = 0; z < 3; z++)
			{
				*(v + 3 * target + z) = *(v_new + 3 * target + z);
			}
			*E += dE;
			si = 1.0;
		}
	}
	return si;
}

double metro_x (double *x, double *p, double *E, double *tabla_V_LJ, double *tabla_V_P, double rc2, double sc2, double dr2, double ds2, double q0, double p0, double L, double beta, double dx, int N)
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
	pos_PBC(x_new, target, L);
	double dE, si;
	dE = variacion_E_x (x, x_new, p, tabla_V_LJ, tabla_V_P, rc2, sc2, dr2, ds2, q0, p0, L, N, target);
	si = aceptacion (x, x_new, E, dE, beta, N, target);
	
	free(x_new);
	return si;
}

double metro_p (double *x, double *p, double *E, double *tabla_V_LJ, double *tabla_V_P, double sc2, double ds2, double q0, double p0, double L, double beta, double dp, int N, double m)
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
	double dE, si;
	dE = variacion_E_p (x, p, p_new, tabla_V_P, sc2, ds2, q0, p0, L, N, target, m);
	si = aceptacion (p, p_new, E, dE, beta, N, target);
	
	free(p_new);
	return si;
}
