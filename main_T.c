#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "tablas.h"
#include "interaccion.h"
#include "metropolis.h"
double deltas (double si_x, double dx, int steps);

int main(int argc,char *argv[])
{
//------Datos de la simulación


	int N = 512; //tiene que ser un cubo y multiplo de 4!!
	double T0 = 4.05, beta = 1.0 / T0;
	double dx, dp; //diferenciales p/ Montecarlo
//	int correlacion = 1, pasos = (int)(0.4 * N * correlacion), pre_termalizacion = 600, termalizacion = 150;
	int correlacion = 20000, pasos = correlacion * 0.1 * 500, pre_termalizacion = 30000000, termalizacion = 500 * 600;
	int l, laps = 20;
	double Temp = T0, dT = -0.04, Tf = 0.05;
	int tfinal =(int)((Tf - Temp) / dT) + 1;
	printf("\n N = %d, correlación = %d, pasos = %d, termalización = %d, pre_termalización = %d, temperaturas = % d\n", N, correlacion, pasos, termalizacion, pre_termalizacion, tfinal);
	double *x, *p;
	x = (double*) malloc(3 * N * sizeof(double));
	p = (double*) malloc(3 * N * sizeof(double));
//--tipo de partículas:-->//(0:N/2) n y (N/2:N) p
//---ctes de los potenciales
	double m = 938.3;		//MeV; masa de p y n (la = p/ ambos)
	double p0 = 2.067 * pow(10, -22); // MeV * s / fm
	double q0 = 6.0; //fm
//---------Tablas de Energía
	int largo_tabla = 500000;
	double *tabla_V_LJ, *tabla_V_P;
	tabla_V_LJ = (double*) malloc(largo_tabla * sizeof(double));
	tabla_V_P = (double*) malloc(largo_tabla * sizeof(double));
	double rc2 = 5.4 * 5.4;
	double sc2 = 10.0;
	double dr2 = tablas (tabla_V_LJ, tabla_V_P, largo_tabla, rc2, sc2, q0, p0);
	double ds2 = dr2 * sc2 / rc2;

	double *E, *Energia;
	E = (double*) malloc(1 * sizeof(double));
	Energia = (double*) malloc(pasos / correlacion * sizeof(double));

	double *Energy, *aceptacion_x, *aceptacion_p;
	Energy = (double*) malloc(tfinal * sizeof(double));
	aceptacion_x = (double*) malloc(tfinal * sizeof(double));
	aceptacion_p = (double*) malloc(tfinal * sizeof(double));

	int t, n;
	double va, total = (double) (tfinal * pasos + laps - 2);
	double si_x, si_p;

	int z;
	double L, drho = -0.01, rho = 0.13 - drho;
	for (z = 0; z < 2; z++)
	{
		rho += drho;
		printf("\nrho ---> %.2lf\t", rho);
		L = cbrt(N / rho);
	//------condiciones iniciales------(arreglo cúbico en las posiciones y gaussiana en los momentos
		T0 = 4.05;
		Temp = 4.05;
		srand(1.0);
		set_pos(x, N, L);
		set_momentos(p, N, T0);
	//---Termalizacion inicial
		*E = hamiltoneano (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, m);
		printf("Energía inicial %lf\n", *E / (double) N);
		dp = 1.2;
		dx = 0.2;
		for (l = 0; l < laps; l++)
		{
			si_x = 0.0;
			si_p = 0.0;
			for (n = 0; n < pre_termalizacion / laps; n++)
			{
			si_x += metro_x (x, p, E, tabla_V_LJ, tabla_V_P, rc2, sc2, dr2, ds2, q0, p0, L, beta, dx, N);
			si_p += metro_p (x, p, E, tabla_V_LJ, tabla_V_P, sc2, ds2, q0, p0, L, beta, dp, N, m);
			}
				va = (double) l * 100.0 / total;
				printf("Progreso %.2lf", va);
				printf("%%\r");

				dx = deltas (si_x, dx, pre_termalizacion / laps);
				dp = deltas (si_p, dp, pre_termalizacion / laps);
		}
	//--------Loop de temperaturas

		for (t = 0; t < tfinal; t++)
		{
			for (l = 0; l < laps; l++)
			{
				si_x = 0.0;
				si_p = 0.0;
				for (n = 0; n < termalizacion / laps; n++)
				{
					si_x += metro_x (x, p, E, tabla_V_LJ, tabla_V_P, rc2, sc2, dr2, ds2, q0, p0, L, beta, dx, N);
					si_p += metro_p (x, p, E, tabla_V_LJ, tabla_V_P, sc2, ds2, q0, p0, L, beta, dp, N, m);
				}
				dx = deltas (si_x, dx, termalizacion / laps);
				dp = deltas (si_p, dp, termalizacion / laps);
			}
			si_x = 0.0;
			si_p = 0.0;
			for (n = 0; n < pasos; n++)
			{
				si_x += metro_x (x, p, E, tabla_V_LJ, tabla_V_P, rc2, sc2, dr2, ds2, q0, p0, L, beta, dx, N);
				si_p += metro_p (x, p, E, tabla_V_LJ, tabla_V_P, sc2, ds2, q0, p0, L, beta, dp, N, m);

				if(n % correlacion == 0)
				{
					*(Energia + n / correlacion) = *E;
				}
				va = (double) (n + t * pasos + laps - 1) * 100.0 / total;

				printf("Progreso %.2lf", va);
				printf("%%\r");
			}
			*(Energy + t) = mean (Energia, 0, pasos / correlacion, 1);
			*(aceptacion_x + t) = si_x / (double) pasos * 100.0;
			*(aceptacion_p + t) = si_p / (double) pasos * 100.0;
			Temp += dT;
			beta = 1.0 / Temp;
		}
		FILE * fp;
		char filename[500];
		sprintf (filename,"/home/pedro/Desktop/Metropolis_Pauli/Data/512F/%d_barrido_rho_%.2lf.txt", N, rho);
		fp = fopen(filename, "w");
		for (t = 0; t < tfinal; t++)
		{
			fprintf(fp, "%lf\t", T0);
			fprintf(fp, "%lf\t", *(Energy + t) / (double) N);
			fprintf(fp, "%lf\t", *(aceptacion_x + t));
			fprintf(fp, "%lf\n", *(aceptacion_p + t));
			T0 += dT;
		}
		fclose(fp);
	}
	free(x);
	free(p);
	free(tabla_V_LJ);
	free(tabla_V_P);
	free(E);
	free(Energia);
	free(Energy);
	free(aceptacion_x);
	free(aceptacion_p);
	return 0;
}

double deltas (double si_x, double dx, int steps)
{	
	if (si_x * 100.0 / (double) steps <= 45.0)
	{
		dx -= dx * 0.1;
	}
	else if (si_x * 100.0 / (double) steps >= 55.0)
	{
		dx += dx * 0.1;
	}
	return dx;
}
#include "general.c"
#include "inicializar.c"
#include "tablas.c"
#include "interaccion.c"
#include "metropolis.c"
