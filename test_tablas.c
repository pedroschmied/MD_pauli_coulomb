#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "tablas.h"

int main()
{
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
//---------
	int i;
	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Final_compu/Datos/test_tablas_V_P.txt");
	fp = fopen(filename, "w");

	for (i = 10000; i < largo_tabla; i++)
	{
		fprintf(fp, "%lf\t", ds2 * (i + 1.0));
		fprintf (fp, "%.45lf\n",  *(tabla_V_P + i));
	}
	fclose(fp);
	double rij2 = 4.0;
	double s2 = 4.5;
	double V_LJ = interpol (tabla_V_LJ, rij2, dr2);
	double V_P = interpol (tabla_V_P, s2, ds2);
	printf("\n Interpolo tablas\n");
	printf("\t Lennard Jones con rij^2 = %lf", rij2);
	printf("----> V_LJ = %lf\n", V_LJ);
	printf("\t Pauli con s^2 = %lf\t", s2);
	printf("-------> V_P = %lf\n", V_P);
	free(tabla_V_LJ);
	free(tabla_V_P);

	return 0;
}

#include "tablas.c"

