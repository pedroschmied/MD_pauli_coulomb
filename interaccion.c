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

int pos_PBC(double *x, int i, double L)
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

double hamiltoneano (double *x, double *p, double *tabla_V_LJ, double *tabla_V_P, double dr2, double ds2, double rc2, double sc2, double L, int N, double q0, double p0, double m)
{
	int i, j;
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));

	double H, rij2, pij2;

	double V_LJ = 0.0, V_P = 0.0, E_cin = 0.0;
	for (i = 0; i < 3 * N; i++)
	{
		E_cin += *(p + i) * *(p + i);
	}
	E_cin = E_cin / (2.0 * m);
	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			delta_x(x, i, x, j, delta_r);
			interaccion_PCB (delta_r, L);
			rij2 = norma2 (delta_r);
//---evalúo L-J
			V_LJ += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);

			delta_x (p, i, p, j, delta_r);
			pij2 = norma2 (delta_r);
//---evalúo Pauli
			if((int)((double)i / (N / 4.0)) == (int)((double)j / (N / 4.0)))
			{
				V_P += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
			}
		}
	}
	
	H = V_LJ + V_P + E_cin;

	free(delta_r);
	return H;
}

double eval_LJ (double *tabla_V_LJ, double rij2, double rc2, double dr2)
{
	double V_LJ = 0.0;
	if (rij2 < rc2)
	{
		V_LJ = interpol (tabla_V_LJ, rij2, dr2);
	}
	return V_LJ;
}

double eval_P (double *tabla_V_P, double rij2, double pij2, double q0, double p0, double sc2, double ds2)
{
	double V_P = 0.0, s2 = pij2 / (p0 * p0) + rij2 / (q0 * q0);	
	if (s2 < sc2)
	{
		V_P = interpol (tabla_V_P, s2, ds2);
	}
	return V_P;
}

double  variacion_E_x (double *x, double *x_new, double *p, double *tabla_V_LJ, double *tabla_V_P, double rc2, double sc2, double dr2, double ds2, double q0, double p0, double L, int N, int target)
{
	int j;
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));
	double H = 0.0, H_new = 0.0, rij2, pij2;

	int b = (int)((double)target / (N / 4.0));
	int a2 = (int)(((double)b + 1.0) * (N / 4.0)), a1 = a2 -(int)(N / 4.0);
//------contribución energética del target vs del taget modificado
	for (j = 0; j < a1; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
//--variación en x
		delta_x(x_new, target, x_new, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H_new += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
	}
	for (j = a1; j < target; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
		delta_x (p, target, p, j, delta_r);
		pij2 = norma2 (delta_r);
		H += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
//--variación en x
		delta_x(x_new, target, x_new, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H_new += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
		H_new += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
	}
	for (j = target + 1; j < a2; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
		delta_x (p, target, p, j, delta_r);
		pij2 = norma2 (delta_r);
		H += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
//--variación en x
		delta_x(x_new, target, x_new, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H_new += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
		H_new += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
	}
	for (j = a2; j < N; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
//--variación en x
		delta_x(x_new, target, x_new, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);
		H_new += eval_LJ (tabla_V_LJ, rij2, rc2, dr2);
	}

	double dE = H_new - H;
	free(delta_r);
	return dE;
}

double  variacion_E_p (double *x, double *p, double *p_new, double *tabla_V_P, double sc2, double ds2, double q0, double p0, double L, int N, int target, double m)
{
	int j;
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));
	double H = 0.0, H_new = 0.0, rij2, pij2;

	int b = (int)((double)target / (N / 4.0));
	int a2 = (int)(((double)b + 1.0) * (N / 4.0)), a1 = a2 -(int)(N / 4.0);
//------contribución energética del target vs del taget modificado
// si no varío x ----> V_LJ_new - V_LJ = 0;
	for (j = a1; j < target; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);

		delta_x (p, target, p, j, delta_r);
		pij2 = norma2 (delta_r);
		H += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
//--variación en x
		delta_x (p_new, target, p_new, j, delta_r);
		pij2 = norma2 (delta_r);
		H_new += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
	}
	for (j = target + 1; j < a2; j++)
	{
		delta_x(x, target, x, j, delta_r);
		interaccion_PCB (delta_r, L);
		rij2 = norma2 (delta_r);

		delta_x (p, target, p, j, delta_r);
		pij2 = norma2 (delta_r);
		H += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
//--variación en x
		delta_x (p_new, target, p_new, j, delta_r);
		pij2 = norma2 (delta_r);
		H_new += eval_P (tabla_V_P, rij2, pij2, q0, p0, sc2, ds2);
	}

	double dE_cin = 0.0;
	for (j = 0; j < 3 * N; j++)
	{
		dE_cin += *(p_new + j) * *(p_new + j) - *(p + j) * *(p + j);
	}
	dE_cin = dE_cin / (2.0 * m);
	double dE = H_new - H + dE_cin;
	free(delta_r);
	return dE;
}
