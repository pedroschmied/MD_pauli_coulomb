#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "tablas.h"
#include "interaccion.h"
#include "metropolis.h"
int ID_random(double *type, int N);
int ID_half(double *type, int N);
int imprimir_ID(double *type, int N);


int main()
{
//------condiciones iniciales------(arreglo cúbico en las posiciones y gaussiana en los momentos=
	int N = 216; //tiene que ser un cubo!!!!
	double rho = 0.16, L = cbrt(N / rho), T0 = 4.0;
	double beta = 1.0 / T0;
	int correlacion = 1000, pasos = 100 * correlacion;
	double *x, *p, *type;
	x = (double*) malloc(3 * N * sizeof(double));
	p = (double*) malloc(3 * N * sizeof(double));
	type = (double*) malloc(N * sizeof(double));
	srand(1.0);
	set_pos(x, N, L);
	set_momentos(p, N, T0);
//--tipo de partículas: 1--->neutrón y 2---->protón con +/- su proyección de spin
	ID_half(type, N);
	imprimir_ID(type, N);
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
	*E = hamiltoneano (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, type);
	printf("\n Energía inicial %lf\n", *E);
//-------Metrópolis
	double aceptacion_x = 0.0, aceptacion_p = 0.0, dx = 0.2, dp = 1.2;
	int t;
	double va;
	for (t = 0; t < pasos; t++)
	{
		aceptacion_x += metro_x (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dx, beta, type);
		aceptacion_p += metro_p (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dp, beta, type);
		if(t % correlacion == 0)
		{
			*(Energia + t / correlacion) = *E;
		}
		va = (double) t * 100.0 / (double) pasos;
		printf("Progreso %.2lf", va);
		printf("%%\r");
	}
	printf("\n Aceptación en x---> %lf", aceptacion_x / (double) pasos * 100.0);
	printf("%%\n");
	printf("\n Aceptación en p---> %lf", aceptacion_p / (double) pasos * 100.0);
	printf("%%\n");
	printf("\n Aceptación---> %lf", (aceptacion_x + aceptacion_p) / (2.0 * (double) pasos) * 100.0);
	printf("%%\n");
	printf("\n Energía final %lf\n", *E);
//-------
	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Final_compu/Datos/prueba2.txt");
	fp = fopen(filename, "w");

	for (t = 0; t < pasos / correlacion; t++)
	{
		fprintf(fp, "%d\t", t);
		fprintf (fp, "%lf\n",  *(Energia + t));
	}
	fclose(fp);

	free(x);
	free(p);
	free(type);
	free(tabla_V_LJ);
	free(tabla_V_P);
	free(E);
	free(Energia);
	return 0;
}
int ID_random(double *type, int N)
{
	double al, bl;
	int i;
	for(i = 0; i < N; i++)
	{
		al = rand()%2 + 1.0; //tira 2s y 1s
		bl = rand()%2 * 2.0 - 1.0; //tira -1s y 1s
		*(type + i) = bl * al;
	}
	return 0;
}

int ID_half(double *type, int N) //1/2 n & p, y 1/2 up & down
{
	int i;
	for(i = 0; i < N; i = i + 2)
	{
		*(type + i) = 2.0;
	}
	for(i = 1; i < N; i = i + 2)
	{
		*(type + i) = 1.0;
	}
	for(i = 0; i < N; i = i + 4)
	{
		*(type + i) = -*(type + i);
	}
	for(i = 1; i < N; i = i + 4)
	{
		*(type + i) = -*(type + i);
	}
	return 0;
}

int imprimir_ID(double *type, int N)
{
	int p_u = 0, p_d = 0, n_u = 0, n_d = 0, i;
	for(i = 0; i < N; i++)
	{
		if(*(type + i) == -1.0)
		{
			n_d += 1;
		}
		else if(*(type + i) == -2.0)
		{
			p_d += 1;
		}

		else if(*(type + i) == 1.0)
		{
			n_u += 1;
		}
		p_u = N - p_d - n_u - n_d;
	}
	printf("\n--------------\n");
	printf("protón:  up %d", p_u);
	printf("   +   ");
	printf("down %d", p_d);
	printf("-----> %d\n", p_d + p_u);
	printf("neutrón: up %d", n_u);
	printf("   +   ");
	printf("down %d", n_d);
	printf("-----> %d\n", n_d + n_u);
	printf("--------------\n");
	return 0;
}
#include "general.c"
#include "inicializar.c"
#include "tablas.c"
#include "interaccion.c"
#include "metropolis.c"
