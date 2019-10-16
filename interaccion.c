#include "interaccion.h"
#include "general.h"
#include "tablas.h"

int interaccion_PCB(double *delta_r, double L)
{
	int k;
	for (k = 0; k < 3; k++)
	{
		if(*(delta_r + k) > L / 2.0)
		{
			*(delta_r + k) -= L;
		}
		else if(*(delta_r + k) <= - L / 2.0)
		{
			*(delta_r + k) += L;
		}
	}
	return 0;
}

int pos_PBC(double *x, int i, int N, double L)
{
	int z;
	for(z = 0; z < 3; z++)
	{
		if (*(x + 3 * i + z) > L)
		{
			*(x + 3 * i + z) = *(x + 3 * i + z) - L;
		}
		else if (*(x + 3 * i + z) < 0)
		{
			*(x + 3 * i + z) = *(x + 3 * i + z) + L;
		}
	}
	return 0;
}

double hamiltoneano (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double *type)
{
	int i, j;
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));
	double *V_LJ;
	V_LJ = (double*) malloc(1 * sizeof(double));
	double *V_P;
	V_P = (double*) malloc(1 * sizeof(double));

	double H, rij2, pij2;

	*V_LJ = 0.0;
	*V_P = 0.0;
	double E_cin = 0.0;
	for (i = 0; i < 3 * N; i++)
	{
		E_cin += *(p + i) * *(p + i);
	}
	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			delta_x(x, i, x, j, delta_r);
			interaccion_PCB (delta_r, L);
			rij2 = norma2 (delta_r);
			delta_x (p, i, p, j, delta_r);
			pij2 = norma2 (delta_r);
			eval_interaccion (V_LJ, V_P, tabla_V_LJ, tabla_V_P, rij2, pij2, p0, q0, rc2, sc2, ds2, dr2, type, i, j);
		}
	}
	H = (*V_LJ + *V_P) / (double) N + E_cin / (double) N;

	free(delta_r);
	free(V_LJ);
	free(V_P);
	return H;
}

double eval_interaccion (double *V_LJ, double *V_P, double *tabla_V_LJ, double *tabla_V_P, double rij2, double pij2, double p0, double q0, double rc2, double sc2, double ds2, double dr2, double *type, int i, int j)
{
	double s2 = pij2 / (p0 * p0) + rij2 / (q0 * q0);
	if (rij2 < rc2)
	{
		*V_LJ += interpol (tabla_V_LJ, rij2, dr2);
	}
	if (s2 < sc2)
	{
		if(*(type + i) / *(type + j) == 1.0)
		{
			*V_P += interpol (tabla_V_P, s2, ds2);
		}
	}
	return 0.0;
}

double  variacion_E  (double *x, double *p, double *x_new, double *p_new, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, int target, double *type)
{
	int i, j;
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));
	double *V_LJ;
	V_LJ = (double*) malloc(1 * sizeof(double));
	double *V_P;
	V_P = (double*) malloc(1 * sizeof(double));

	double H, H_new, rij2, pij2;
//------contribución energética del target
	*V_LJ = 0.0;
	*V_P = 0.0;
	for (j = 0; j < target; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		delta_x (p, target, p, j, delta_r);
		pij2 = norma2 (delta_r);
		eval_interaccion (V_LJ, V_P, tabla_V_LJ, tabla_V_P, rij2, pij2, p0, q0, rc2, sc2, ds2, dr2, type, target, j);
	}
	for (j = target + 1; j < N; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		delta_x (p, target, p, j, delta_r);
		pij2 = norma2 (delta_r);
		eval_interaccion (V_LJ, V_P, tabla_V_LJ, tabla_V_P, rij2, pij2, p0, q0, rc2, sc2, ds2, dr2, type, target, j);
	}
	H = (*V_LJ + *V_P) / (double) N;
//------contribución energética del target modificado
	*V_LJ = 0.0;
	*V_P = 0.0;
	for (j = 0; j < target; j++)
	{
		delta_x(x_new, target, x_new, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		delta_x (p_new, target, p_new, j, delta_r);
		pij2 = norma2 (delta_r);
		eval_interaccion (V_LJ, V_P, tabla_V_LJ, tabla_V_P, rij2, pij2, p0, q0, rc2, sc2, ds2, dr2, type, target, j);
	}

	for (j = target + 1; j < N; j++)
	{
		delta_x(x_new, target, x_new, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		delta_x (p_new, target, p_new, j, delta_r);
		pij2 = norma2 (delta_r);
		eval_interaccion (V_LJ, V_P, tabla_V_LJ, tabla_V_P, rij2, pij2, p0, q0, rc2, sc2, ds2, dr2, type, target, j);
	}

	H_new = (*V_LJ + *V_P) / (double) N;
//---variación en la energía cinética
	double dE_cin = 0.0;
	for (i = 0; i < 3 * N; i++)
	{
		dE_cin += (*(p_new + i) * *(p_new + i) - *(p + i) * *(p + i)) / (double) N;
	}
	double dE = H_new - H + dE_cin;

	free(delta_r);
	free(V_LJ);
	free(V_P);
	return dE;
}
