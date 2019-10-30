#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "tablas.h"
#include "interaccion.h"
#include "metropolis.h"

int main(int argc,char *argv[])
{
//------Datos de la simulación
	int N = 512; //tiene que ser un cubo y multiplo de 4!!
	double rho;
	sscanf(argv[1],"%lf", &rho);
	printf("rho ---> %.2lf\n", rho);
	double L = cbrt(N / rho), T0 = 4.0;
	double beta = 1.0 / T0;
	double dx, dp; //diferenciales p/ Montecarlo
	int correlacion = 1 ,pasos = 10 * N * correlacion + 1000000, termalizacion = 0;

	double *x, *p;
	x = (double*) malloc(3 * N * sizeof(double));
	p = (double*) malloc(3 * N * sizeof(double));
//------condiciones iniciales------(arreglo cúbico en las posiciones y gaussiana en los momentos
	srand(1.0);
	set_pos(x, N, L);
	set_momentos(p, N, T0);
//--tipo de partículas:-->//(0:N/2) n y (N/2:N) p
//---ctes de los potenciales
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
//-------Energia inicial
	double *E;
	E = (double*) malloc(1 * sizeof(double));
	double *Energia;
	Energia = (double*) malloc(pasos / correlacion * sizeof(double));
	*E = hamiltoneano (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0);
	printf("\n Energía inicial %lf\n", *E / (double) N);
//-------Metrópolis
	int n;
	double va;
	double total = (double) (pasos + termalizacion);
	double ax = 0.146624, bx = -0.013726, ap = 1.02777, bp =  0.0117; //parametros del ajuste para calcular los dx y dp
	dp = sqrt(T0) * ap + bp;
	dx = sqrt(T0) * ax + bx;
	double dep_rho = 0.13 / rho;
	dx = dx * dep_rho;
//---Termalizacion inicial
	double si_x, si_p;
	si_x = 0.0;
	si_p = 0.0;
	for (n = 0; n < pasos; n++)
	{
		si_x += metro_x (x, p, E, tabla_V_LJ, tabla_V_P, rc2, sc2, dr2, ds2, q0, p0, L, beta, dx, N);
		si_p += metro_p (x, p, E, tabla_V_LJ, tabla_V_P, sc2, ds2, q0, p0, L, beta, dp, N);

		if(n % correlacion == 0)
		{
			*(Energia + n / correlacion) = *E;
		}
		va = (double) (n + termalizacion) * 100.0 / total;
		printf("Progreso %.2lf", va);
		printf("%%\r");
	}
	si_x = si_x / (double) pasos * 100.0;
	si_p = si_p / (double) pasos * 100.0;
	printf("\n aceptación en x---->%.2lf\n", si_x);
	printf("\n aceptación en p---->%.2lf\n", si_p);
	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Metropolis_Pauli/Datos/corrida_T_%.2lf_rho_%.2lf.txt", T0, rho);
	fp = fopen(filename, "w");

	for (n = 0; n < pasos; n++)
	{		
		fprintf(fp, "%d\t", n);
		fprintf (fp, "%lf\n",  *(Energia + n) / (double) N);
		va = (double) (total - pasos + n) * 100.0 / total;
		printf("Progreso %.2lf", va);
		printf("%%\r");

	}
	fclose(fp);

	free(x);
	free(p);
	free(tabla_V_LJ);
	free(tabla_V_P);
	free(E);
	free(Energia);
	return 0;
}
#include "general.c"
#include "inicializar.c"
#include "tablas.c"
#include "interaccion.c"
#include "metropolis.c"
