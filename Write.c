/**
 * File              : Write.c
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 29.04.2019
 * Last Modified Date: 10.05.2019
 */
#include "Write.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

void logger(Data* data, unsigned long int n)
{
	char stru[40];
	char strv[40];
	char strt[40];
	char strca[40];
	char strcb[40];
	char strd[40];
	char strw[40];

	unsigned long int nx1 = data->nx1;
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;

	double *u = data->u;
	double *v = data->v;
	double *T = data->T;
	double *Ca = data->Ca;
	double *Cb = data->Cb;

	sprintf(stru, "results/U-%lu.txt", n);
	sprintf(strv, "results/V-%lu.txt", n);
	sprintf(strt, "results/T-%lu.txt", n);
	sprintf(strd, "results/D-%lu.txt", n);
	sprintf(strw, "results/W-%lu.txt", n);
	sprintf(strca, "results/Ca-%lu.txt", n);
	sprintf(strcb, "results/Cb-%lu.txt", n);

	FILE *fu = fopen(stru, "w");
	FILE *fv = fopen(strv, "w");
	FILE *ft = fopen(strt, "w");
	FILE *fd = fopen(strd, "w");
	FILE *fw = fopen(strw, "w");
	FILE *fca = fopen(strca, "w");
	FILE *fcb = fopen(strcb, "w");
	FILE *track = NULL;
   	if (n ==0)
	{
		track = fopen("results/out.log", "w");
		fprintf(track, "div\t");
		fprintf(track, "Tmean\t");
		fprintf(track, "Qub\t");
		fprintf(track, "H\t");
		fprintf(track, "Qb\t");
		fprintf(track, "Reh\t");
		fprintf(track, "Rew\t");
		fprintf(track, "Pet\t");
		fprintf(track, "PeCa\t");
		fprintf(track, "PeCb\n");
	}
	else
		track = fopen("results/out.log", "a");

	double gT, gCa, gCb, uv, w, div;
	double uvmax = -DBL_MAX; 
	double wmax = -DBL_MAX; 
	double gTmax = -DBL_MAX;
	double gCamax = -DBL_MAX;
	double gCbmax = -DBL_MAX;
	
	double sdiv = 0.;
	double sT = 0.;
	double Qb = 0.;
	double Qub = 0.;
	double H = 0.;
	double hinv = data->hinv;
	double h = 1./hinv;
	unsigned long int i, j;
	for (j = 0; j < ny; ++j)
	{
		for (i = 0; i < nx; ++i)
		{
			fprintf(fu, "%f\t", u[i*ny + j]);
			fprintf(fv, "%f\t", v[i*ny + j]);
			uv = fabs(u[i*ny + j]) + fabs(v[i*ny + j]);
			uvmax = (uv > uvmax) ? uv : uvmax;

			fprintf(ft, "%f\t", T[i*ny + j]);
			fprintf(fca, "%f\t", Ca[i*ny + j]);
			fprintf(fcb, "%f\t", Cb[i*ny + j]);
			if ((i > 0 && j > 0) && (i < nx-1 && j < ny-1))
			{ 
				gT = fmax(T[i*ny + j+1] - T[i*ny + j-1], T[(i+1)*ny + j] - T[(i-1)*ny + j]);
				gCa = fmax(Ca[i*ny + j+1] - Ca[i*ny + j-1], Ca[(i+1)*ny + j] - Ca[(i-1)*ny + j]);
				gCb = fmax(Cb[i*ny + j+1] - Cb[i*ny + j-1], Cb[(i+1)*ny + j] - Cb[(i-1)*ny + j]);
				gTmax = (gT > gTmax) ? gT : gTmax;
				gCamax = (gCa > gCamax) ? gCa : gCamax;
				gCbmax = (gCb > gCbmax) ? gCb : gCbmax;

				sT += T[i*ny + j] - 1.;

				div = hinv*(- u[(i-1)*ny + j] + u[i*ny + j] - v[i*ny + j-1] + v[i*ny + j]);
				sdiv += div;
				fprintf(fd, "%f\t", div);

				w = (v[(i+1)*ny + j] - v[i*ny + j] - u[i*ny + j+1] + u[i*ny + j]);
				fprintf(fw, "%f\t", hinv*w);
				wmax = (w > wmax) ? w : wmax;
			}
		}	
		fprintf(fu, "\n");
		fprintf(fv, "\n");
		fprintf(ft, "\n");
		fprintf(fd, "\n");
		fprintf(fw, "\n");
		fprintf(fca, "\n");
		fprintf(fcb, "\n");
	}

	for (i = nx1 + 1; i < 2*nx1; ++i)
	{
		Qub += .5*v[i*ny + ny-2]*(Cb[i*ny + ny-1] + Cb[i*ny + ny-2])/hinv;
		H += .5*v[i*ny + ny-2]*(T[i*ny + ny-1] + T[i*ny + ny-2] - 2.)/hinv;
	}
	for (i = 1; i < nx-1; ++i)
		Qb += .5*k(data, .5*(T[i*ny] + T[i*ny + 1]))*(Ca[i*ny] + Ca[i*ny + 1])/hinv;

	fprintf(track, "%1.5e\t", sdiv);
	fprintf(track, "%1.5e\t", hinv*hinv*sT/18.);
	fprintf(track, "%1.5e\t", Qub);
	fprintf(track, "%1.5e\t", H);
	fprintf(track, "%1.5e\t", Qb);
	fprintf(track, "%1.5e\t", h*uvmax/data->nu);
	fprintf(track, "%1.5e\t", h*wmax/data->nu);
	fprintf(track, "%1.5e\t", .5*h*gTmax/data->alpha);
	fprintf(track, "%1.5e\t", .5*h*gCamax/data->D);
	fprintf(track, "%1.5e\n", .5*h*gCbmax/data->D);

	fclose(fu);
	fclose(fv);
	fclose(ft);
	fclose(fd);
	fclose(fw);
	fclose(fca);
	fclose(fcb);
	fclose(track);
}
