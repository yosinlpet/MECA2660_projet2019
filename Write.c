/**
 * File              : Write.c
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 29.04.2019
 * Last Modified Date: 09.05.2019
 */
#include "Write.h"
#include <stdio.h>
#include <string.h>

void logger(Data* data, unsigned long int n)
{
	char stru[20];
	char strv[20];
	char strp[20];
	char strd[20];
	char strw[20];
	char strf[20];

	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;

	double *u = data->u;
	double *v = data->v;
	double *p = data->p;
	double *phi = data->phi;

	sprintf(stru, "results/U-%lu.txt", n);
	sprintf(strv, "results/V-%lu.txt", n);
	sprintf(strp, "results/P-%lu.txt", n);
	sprintf(strd, "results/D-%lu.txt", n);
	sprintf(strw, "results/W-%lu.txt", n);
	sprintf(strf, "results/F-%lu.txt", n);

	FILE *fu = fopen(stru, "w");
	FILE *fv = fopen(strv, "w");
	FILE *fp = fopen(strp, "w");
	FILE *fd = fopen(strd, "w");
	FILE *fw = fopen(strw, "w");
	FILE *ff = fopen(strf, "w");

	double wmax = 0.; 
	double sdiv = 0.;
	double hinv = data->hinv;
	double h = 1./hinv;
	double w, div;
	unsigned long int i, j;
	for (j = 0; j < ny; ++j)
	{
		for (i = 0; i < nx; ++i)
		{
			fprintf(fu, "%f\t", u[i*ny + j]);
			fprintf(fv, "%f\t", v[i*ny + j]);
			fprintf(fp, "%f\t", p[i*ny + j]);
			fprintf(ff, "%f\t", phi[i*ny + j]);
			if ((i > 0 && j > 0) && (i < nx-1 && j < ny-1))
			{ 
				div = hinv*(- u[(i-1)*ny + j] + u[i*ny + j] - v[i*ny + j-1] + v[i*ny + j]);
				sdiv += div;
				fprintf(fd, "%f\t", div);

				w = (v[(i+1)*ny + j] - v[i*ny + j] - u[i*ny + j+1] + u[i*ny + j]);
				fprintf(fw, "%f\t", w);
				wmax = (w > wmax) ? w : wmax;
			}
		}	
		fprintf(fu, "\n");
		fprintf(fv, "\n");
		fprintf(fp, "\n");
		fprintf(fd, "\n");
		fprintf(fw, "\n");
		fprintf(ff, "\n");
	}
	printf("Rew_MAX = %f\n", h*wmax/data->nu);
	printf("Sum Divergence = %1.3e\n", sdiv);

	fclose(fu);
	fclose(fv);
	fclose(fp);
	fclose(fd);
	fclose(fw);
	fclose(ff);
}
