#include "tablas.h"

double tablas (double *tabla_V_LJ, double *tabla_V_P, int largo_tabla, double rc2, double sc2, double q0, double p0)
{
//________ctes______
	double V_P0 = 34.32, V_LJ0 = 18.263; //MeV
	double r1 = 1.7456, r2 = 1.7324, d = 3.35; //fm
	double p1 = 6.2, p2 = 3.0;
	double alfa = 1.2; //fm^(-1)
	double h_barra = 6.582119624 * pow(10, -22); // MeV * s
//__________________
	double H = (h_barra / p0) * 1.0 / q0, cte_pauli = V_P0 * H * H * H;


//---------Lennard-Jones----
	double dr2 = rc2 / ((double)largo_tabla + 1.0), ds2 = sc2 / ((double)largo_tabla + 1.0); //así tengo "largo_tabla" ptos y el r02 = dr2 para simplificar cuentas.
	double V_LJ_rc2 = V_LJ0 * (pow(r1 * r1 / rc2, p1 / 2.0) - pow(r2 * r2 / rc2, p2 / 2.0)) / (1.0 + exp(alfa * (sqrt(rc2) - d)));
	double rij2 = 0.0, s2 = 0.0;
	int l;
	for (l = 0; l < largo_tabla; l++)
	{
		rij2 += dr2;
		*(tabla_V_LJ + l) = V_LJ0 * (pow(r1 * r1 / rij2, p1 / 2.0) - pow(r2 * r2 / rij2, p2 / 2.0)) / (1.0 + exp(alfa * (sqrt(rij2) - d))) - V_LJ_rc2;

		s2 += ds2;
		*(tabla_V_P + l) = cte_pauli * exp(-s2 / 2.0) - cte_pauli * exp(-sc2 / 2.0);
	}

//-------------Pauli// s2 = p^2 / (p0^2) + r^2 / (q0^2)
	// ds = dr2 * sc2 / rf2 = dr2 * 10.0 / rc^2;
	return dr2;
}

double interpol (double *tabla, double rij, double dr)
{
	int k = (int)((rij - dr) / dr); // r0 = dr
	double a;
	if(k < 0)
	{
		k = 0;
	}
	a = (*(tabla + k + 1) - *(tabla + k)) * (rij - (double) k * dr) / dr + *(tabla + k);
	return a;
}
