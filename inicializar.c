#include "general.h"
#include "inicializar.h"
//	distribuyo la posicón de las partículas en una caja 3D, elegí (x0, y0, z0) = (0.5, 0.5, 0.5) solo para evitar el origen xq Lennard Jones diverge allí.
double set_pos(double *x, int N, double L)
{
	int n = cbrt(N);
	int i = 0, x1, x2, x3;
	double dl = L /(double) n;
	for(x1 = 0; x1 < n; x1++)
	{
		for(x2 = 0; x2 < n; x2++)
		{
			for(x3 = 0; x3 < n; x3++)
			{
				*(x + 3 * i + 0) = dl * ((double) x1 + 0.5);
				*(x + 3 * i + 1) = dl * ((double) x2 + 0.5);
				*(x + 3 * i + 2) = dl * ((double) x3 + 0.5);
				i++;
			}
		}
	}
	return dl;
}
// distribución gaussiana en la velocidad y le resto la Veloc. del centro de masa para que no haya un flujo de partículas
double set_momentos(double *v, int N, double T)
{
	double sigma = sqrt(T);// m = 1
	int i, k;
	srand(100.0);
	for(i = 0; i < 3 * N; i++)
	{
		*(v + i) = gaussiana(0.0, sigma);
	}

	double *vcm;
	vcm = (double*) malloc(3 * sizeof(double));
	for(k = 0; k < 3; k++)
	{
		*(vcm + k) = 0.0;
	}

	for(i = 0; i < N; i++)
	{
		for(k = 0; k < 3; k++)
		{
			*(vcm + k) += *(v + 3 * i + k) /(double) N;
		}
	}
	for(i = 0; i < N; i++)
	{
		for(k = 0; k < 3; k++)
		{
			*(v + 3 * i + k) -= *(vcm + k);
		}
	}
	free(vcm);
	double E_cin = 0.0;
	for(i = 0; i < 3 * N; i++)
	{
		E_cin += *(v + i) * *(v + i) / 2.0;
	}
	return E_cin;
}
