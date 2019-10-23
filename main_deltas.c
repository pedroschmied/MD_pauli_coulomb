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
//------Datos de la simulación
	int N = 512; //tiene que ser un cubo!!!!
	double rho = 0.14, L = cbrt(N / rho), T0 = 4.05;
	double beta = 1.0 / T0;
	double dx, dp; //diferenciales p/ Montecarlo
//	int correlacion = 500, pasos = 1000 * correlacion, pre_termalizacion = 30000000, termalizacion = 50000;20000000
	int l = 50, correlacion = 1, pasos = l * 10 * N * correlacion, pre_termalizacion = 20000000, termalizacion = 5;

	double Temp = T0, dT = -0.025, Tf = 0.05 + dT;
	int tfinal =(int)((Tf - Temp) / dT);
	double *x, *p, *type;
	x = (double*) malloc(3 * N * sizeof(double));
	p = (double*) malloc(3 * N * sizeof(double));
	type = (double*) malloc(N * sizeof(double));
//------condiciones iniciales------(arreglo cúbico en las posiciones y gaussiana en los momentos
	srand(1.0);
	set_pos(x, N, L);
	set_momentos(p, N, T0);
//--tipo de partículas: 1--->neutrón y 2---->protón con +/- su proyección de spin
	ID_half(type, N);
//	imprimir_ID(type, N);
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
	int t, n;
	double va;
	double total = (double) (tfinal + pre_termalizacion);
	dp = 2.0;
	dx = 0.5;
//---Termalizacion inicial
	for (n = 0; n < pre_termalizacion; n++)
	{
		metro_x (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dx, beta, type);
		metro_p (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dp, beta, type);
		va = (double) n * 100.0 / total;
		printf("Progreso %.2lf", va);
		printf("%%\r");
	}
//--------Loop de temperaturas
	double *Energy;
	Energy = (double*) malloc(tfinal * sizeof(double));
	double *aceptacion_x, *aceptacion_p;
	aceptacion_x = (double*) malloc(tfinal * sizeof(double));
	aceptacion_p = (double*) malloc(tfinal * sizeof(double));
	double si_x, si_p;
	double *Dx, *Dp;
	Dx = (double*) malloc(tfinal * sizeof(double));
	Dp = (double*) malloc(tfinal * sizeof(double));
	int h;
	for (t = 0; t < tfinal; t++)
	{
		for (n = 0; n < termalizacion; n++)
		{
			metro_x (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dx, beta, type);
			metro_p (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dp, beta, type);
		}
		for(h = 0; h < l - 1; h++)
		{
			si_x = 0.0;
			si_p = 0.0;
			for (n = 0; n < pasos / l; n++)
			{
				si_x += metro_x (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dx, beta, type);
				si_p += metro_p (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dp, beta, type);
				*(Energia + (n + h * pasos / l) / correlacion) = *E;			

			}	

			if(si_x / (double)(pasos / l) <= 0.45)
			{
				dx -= dx / 10.0; //le bajo un 10%
			}
			if(si_p / (double)(pasos / l) <= 0.45)
			{
				dp -= dp / 10.0; //le bajo un 10%
			}
		}
		si_x = 0.0;
		si_p = 0.0;
		for (n = 0; n < pasos / l; n++)
		{
			si_x += metro_x (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dx, beta, type);
			si_p += metro_p (x, p, tabla_V_LJ, tabla_V_P, dr2, ds2, rc2, sc2, L, N, q0, p0, E, dp, beta, type);
			*(Energia + (n + pasos * (l - 1) / l) / correlacion) = *E;			
		}
		*(Energy + t) = mean (Energia, 0, pasos / correlacion, 1);
		*(aceptacion_x + t) = si_x / (double) (pasos / l) * 100.0;
		*(aceptacion_p + t) = si_p / (double) (pasos / l) * 100.0;
		*(Dx + t) = dx;
		*(Dp + t) = dp;
		
		Temp += dT;
		beta = 1.0 / Temp;
		va = (double) (t + pre_termalizacion) * 100.0 / total;
		printf("Progreso %.2lf", va);
		printf("%%\r");
	}
	printf("\nGuardo datos\n");
	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Final_compu/Datos/barrido_deltas_rho_%.2lf.txt", rho);
	fp = fopen(filename, "w");

	for (t = 0; t < tfinal; t++)
	{		
		fprintf(fp, "%lf\t", T0);
		fprintf (fp, "%lf\t",  *(Energy + t));
		fprintf(fp, "%lf\t",  *(aceptacion_x + t));
		fprintf(fp, "%lf\t",  *(aceptacion_p + t));
		fprintf(fp, "%lf\t",  *(Dx + t));
		fprintf(fp, "%lf\n",  *(Dp + t));
		T0 += dT;
	}
	fclose(fp);

	free(x);
	free(p);
	free(type);
	free(tabla_V_LJ);
	free(tabla_V_P);
	free(E);
	free(Energia);
	free(Energy);
	free(aceptacion_x);
	free(aceptacion_p);
	free(Dx);
	free(Dp);
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
