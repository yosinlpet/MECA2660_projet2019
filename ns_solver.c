/** * File              : ns_solver.c
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 25.04.2019
 * Last Modified Date: 10.05.2019
 */

#include <stdlib.h>
#include <math.h>
#include "ns_solver.h"
#include "poisson.h"

void step_one(Data* data) {
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;
	double hinv = data->hinv;
	double dt = data->dt;
	double nu = data->nu;
	double beta = data->beta;

	double g = data->g;
	double *T = data->T;

	double *p = data->p;
	double *u = data->u;
	double *v = data->v;
	double *old_u = data->old_u;
	double *old_v = data->old_v;

	double *ucopy = malloc(nx*ny*sizeof(double));
	double *vcopy = malloc(nx*ny*sizeof(double));

	unsigned long int i, j;
	for (i = 0; i < nx; ++i)
	{
		for (j = 0; j < ny; ++j)
		{
			ucopy[i*ny + j] = u[i*ny + j];
			vcopy[i*ny + j] = v[i*ny + j];
		}
	}

	for (i = 1; i < nx-2; ++i)
	{
		for (j = 1; j < ny-1; ++j)
		{
			//Advection
			u[i*ny + j] -= .375*hinv*dt*((ucopy[(i+1)*ny + j] + ucopy[i*ny + j])*(ucopy[(i+1)*ny + j] - ucopy[i*ny + j])
							+ (ucopy[i*ny + j] + ucopy[(i-1)*ny + j])*(ucopy[i*ny + j] - ucopy[(i-1)*ny + j])
							+ (vcopy[(i+1)*ny + j] + vcopy[i*ny + j])*(ucopy[i*ny + j+1] - ucopy[i*ny + j])
							+ (vcopy[(i+1)*ny + j-1] + vcopy[i*ny + j-1])*(ucopy[i*ny + j] - ucopy[i*ny + j-1]));

			u[i*ny + j]  += .125*hinv*dt*((old_u[(i+1)*ny + j] + old_u[i*ny + j])*(old_u[(i+1)*ny + j] - old_u[i*ny + j])
							+ (old_u[i*ny + j] + old_u[(i-1)*ny + j])*(old_u[i*ny + j] - old_u[(i-1)*ny + j])
							+ (old_v[(i+1)*ny + j] + old_v[i*ny + j])*(old_u[i*ny + j+1] - old_u[i*ny + j])
							+ (old_v[(i+1)*ny + j-1] + old_v[i*ny + j-1])*(old_u[i*ny + j] - old_u[i*ny + j-1]));

			//Pressure
			u[i*ny + j] -= hinv*dt*(p[(i+1)*ny + j] - p[i*ny + j]);

			//Diffusion
			u[i*ny + j] += nu*dt*hinv*hinv*(ucopy[(i+1)*ny + j] - 4.*ucopy[i*ny + j] + ucopy[(i-1)*ny + j]
										  + ucopy[i*ny + j+1] + ucopy[i*ny + j-1]);
		}
	}

	for (i = 1; i < nx-1; ++i)
	{
		for (j = 1; j < ny-2; ++j)
		{

			//Advection
			v[i*ny + j] -= .375*hinv*dt*((ucopy[i*ny + j+1] + ucopy[i*ny + j])*(vcopy[(i+1)*ny + j] - vcopy[i*ny + j])
							+ (ucopy[(i-1)*ny + j+1] + ucopy[(i-1)*ny + j])*(vcopy[i*ny + j] - vcopy[(i-1)*ny + j])
							+ (vcopy[i*ny + j+1] + vcopy[i*ny + j])*(vcopy[i*ny + j+1] - vcopy[i*ny + j])
							+ (vcopy[i*ny + j] + vcopy[i*ny + j-1])*(vcopy[i*ny + j] - vcopy[i*ny + j-1]));

			v[i*ny + j] += .125*hinv*dt*((old_u[i*ny + j+1] + old_u[i*ny + j])*(old_v[(i+1)*ny + j] - old_v[i*ny + j])
							+ (old_u[(i-1)*ny + j+1] + old_u[(i-1)*ny + j])*(old_v[i*ny + j] - old_v[(i-1)*ny + j])
							+ (old_v[i*ny + j+1] + old_v[i*ny + j])*(old_v[i*ny + j+1] - old_v[i*ny + j])
							+ (old_v[i*ny + j] + old_v[i*ny + j-1])*(old_v[i*ny + j] - old_v[i*ny + j-1]));

			//Pressure
			v[i*ny + j] -= hinv*dt*(p[i*ny + j+1] - p[i*ny + j]);

			//Diffusion
			v[i*ny + j] += nu*dt*hinv*hinv*(vcopy[i*ny + j+1] - 4.*vcopy[i*ny + j] + vcopy[i*ny + j-1]
										  + vcopy[(i+1)*ny + j] + vcopy[(i-1)*ny + j]);

			//Boussinesq
			v[i*ny + j] -= beta*dt*(T[i*ny + j] - 1.)*g;
		}
	}

	for (i = 0; i < nx; ++i)
	{
		for (j = 0; j < ny; ++j)
		{
			old_u[i*ny + j] = ucopy[i*ny + j];
			old_v[i*ny + j] = vcopy[i*ny + j];
		}
	}

	free(ucopy);
	free(vcopy);
}

void reprojection(Data* data)
{
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;
	double hinv = data->hinv;
	double dt = data->dt;

	double *u = data->u;
	double *v = data->v;
	double *phi = data->phi;
	double *p = data->p;

	unsigned long int i, j;
	for (i = 1; i < nx-2; ++i)
		for (j = 1; j < ny-1; ++j)
			u[i*ny + j] -= hinv*dt*(phi[(i+1)*ny + j] - phi[i*ny + j]);
	for (i = 1; i < nx-1; ++i)
		for (j = 1; j < ny-2; ++j)
			v[i*ny + j] -= hinv*dt*(phi[i*ny + j+1] - phi[i*ny + j]);
	for (i = 1; i < nx-1; ++i)
		for (j = 1; j < ny-1; ++j)
			p[i*ny + j] += phi[i*ny + j];
}

void impose_boundary(Data* data)
{
	unsigned long int nx1 = data->nx1;
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;

	double alpha = data->alpha;
	double D = data->D;
	double DH = data->DH;

	double dt = data->dt;
	double hinv = data->hinv;

	double *v = data->v;
	double *T = data->T;
	double *Ca = data->Ca;
	double *Cb = data->Cb;

	unsigned long int i, j;
	for (j = 1; j < ny-1; ++j) 
	{
		//LEFT
		if (j <= nx1/4 || j > 5*nx1/4)
		{
			T[j] = T[ny + j];
			Ca[j] = Ca[ny + j];
			Cb[j] = Cb[ny + j];
		}

		//RIGHT
		T[(nx-1)*ny + j] = T[(nx-2)*ny + j];
		Ca[(nx-1)*ny + j] = Ca[(nx-2)*ny + j];
		Cb[(nx-1)*ny + j] = Cb[(nx-2)*ny + j];
	}

	double stuff = .5/(D*hinv);
	for (i = 1; i < nx-1; ++i)
	{
		//BOT
		Ca[i*ny] = (1.-stuff*k(data, T[i*ny]))*Ca[i*ny + 1]/(1. + stuff*k(data, T[i*ny]));
		Cb[i*ny] = Cb[i*ny + 1] + k(data, T[i*ny])*stuff*(Ca[i*ny] + Ca[i*ny + 1]);
		T[i*ny] = T[i*ny] + .5*DH*k(data, T[i*ny])*(Ca[i*ny] + Ca[i*ny + 1])/(alpha*hinv);

		//TOP
		T[i*ny + ny-1] = T[i*ny + ny-2];
		Ca[i*ny + ny-1] = Ca[i*ny + ny-2];
		Cb[i*ny + ny-1] = Cb[i*ny + ny-2];
	}

	for (i = 1 + nx1; i <= 2*nx1; ++i) //Don't forget the useless ghost points.
	{
		//TOP
		v[i*ny + ny-2] = v[i*ny + ny-2] - dt*hinv*(v[i*ny + ny-2] - v[i*ny + ny-3]);
		T[i*ny + ny-1] = T[i*ny + ny-1] - dt*hinv*(T[i*ny + ny-1] - T[i*ny + ny-2]);
		Ca[i*ny + ny-1] = Ca[i*ny + ny-1] - dt*hinv*(Ca[i*ny + ny-1] - Ca[i*ny + ny-2]);
		Cb[i*ny + ny-1] = Cb[i*ny + ny-1] - dt*hinv*(Cb[i*ny + ny-1] - Cb[i*ny + ny-2]);
	}
}

double k(Data *data, double T)
{
	return data->k0*exp(-data->ER/T);
}

void ghost_points(Data *data, double vgl[])
{
	unsigned long int nx1 = data->nx1;
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;

	double *u = data->u;
	double *v = data->v;

	unsigned long int i, j;
	for (i = 1; i < nx-2; ++i) //Don't forget the useless ghost points.
	{
		//Bottom
		u[i*ny] = -.2*(u[i*ny + 3] - 5*u[i*ny + 2] + 15*u[i*ny + 1]);

		//Top
		if (i > nx1 && i < 2*nx1)
			u[i*ny + ny-1] = u[i*ny + ny-2]; //Natural outflow
		else
			u[i*ny + ny-1] = -.2*(u[i*ny + ny-4] - 5*u[i*ny + ny-3] + 15*u[i*ny + ny-2]);
	}

	for (j = 1; j < ny-2; ++j)
	{
		//LEFT
		v[j] = -.2*(v[3*ny + j] - 5*v[2*ny + j] + 15*v[ny + j] - 16*vgl[j]); 
		
		//RIGHT
		v[(nx-1)*ny + j] = -.2*(v[(nx-4)*ny + j] - 5*v[(nx-3)*ny + j] + 15*v[(nx-2)*ny + j]);
	}
}

void flow_conservation(Data *data)
{
	unsigned long int nx1 = data->nx1;
	unsigned long int ny = data->ny;

	double hinv = data->hinv;
	double *u = data->u;
	double *v = data->v;

	double Qu = 0.;
	double Qv = 0.;
	unsigned long int j;
	for (j = nx1/4 + 1; j <= 5*nx1/4 ; ++j)
	{
		Qu += u[j];
		Qv += v[(j+3*nx1/4)*ny + ny-2];
	}

	double dv = (Qu - Qv) / hinv;
	for (j = nx1 + 1; j <= 2*nx1 ; ++j)
		v[j*ny + ny-2] += dv;
}

void temp(Data *data)
{
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;

	double alpha = data->alpha;
	double dt = data->dt;
	double hinv = data->hinv;

	double *old_u = data->old_u;
	double *old_v = data->old_v;
	double *old_T = data->old_T;
	double *u = data->u;
	double *v = data->v;
	double *T = data->T;
	double *Tcopy = malloc(nx*ny*sizeof(double));

	unsigned long int i, j;
	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
			Tcopy[i*ny + j] = T[i*ny + j];

	for (i = 1; i < nx-1; ++i)
		for (j = 1; j < ny-1; ++j)
		{
			T[i*ny + j] -= .375*dt*hinv*((u[i*ny + j] + u[(i-1)*ny + j])*(Tcopy[(i+1)*ny + j] - Tcopy[(i-1)*ny + j]) +
										(v[i*ny + j] + v[i*ny + j-1])*(Tcopy[i*ny + j+1] - Tcopy[i*ny + j-1]));

			T[i*ny + j] += .125*dt*hinv*((old_u[i*ny + j] + old_u[(i-1)*ny + j])*(old_T[(i+1)*ny + j] - old_T[(i-1)*ny + j]) +
										(old_v[i*ny + j] + old_v[i*ny + j-1])*(old_T[i*ny + j+1] - old_T[i*ny + j-1]));

			T[i*ny + j] +=  alpha*dt*hinv*hinv*(Tcopy[i*ny + j+1] - 4.*Tcopy[i*ny + j] + Tcopy[i*ny + j-1]
										  + Tcopy[(i+1)*ny + j] + Tcopy[(i-1)*ny + j]);
		}

	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
			old_T[i*ny + j] = Tcopy[i*ny + j];

	free(Tcopy);
}

void concentration(Data *data, char type)
{
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;

	double D = data->D;
	double dt = data->dt;
	double hinv = data->hinv;

	double *old_u = data->old_u;
	double *old_v = data->old_v;
	double *u = data->u;
	double *v = data->v;

	double *C, *old_C;
	if (type == 'a')
	{
		C = data->Ca;
		old_C = data->old_Ca;
	}
	else
	{
		C = data->Cb;
		old_C = data->old_Cb;
	}

	double *Ccopy = malloc(nx*ny*sizeof(double));

	unsigned long int i, j;
	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
			Ccopy[i*ny + j] = C[i*ny + j];

	for (i = 1; i < nx-1; ++i)
		for (j = 1; j < ny-1; ++j)
		{
			C[i*ny + j] -= .375*dt*hinv*((u[i*ny + j] + u[(i-1)*ny + j])*(Ccopy[(i+1)*ny + j] - Ccopy[(i-1)*ny + j]) +
										(v[i*ny + j] + v[i*ny + j-1])*(Ccopy[i*ny + j+1] - Ccopy[i*ny + j-1]));

			C[i*ny + j] += .125*dt*hinv*((old_u[i*ny + j] + old_u[(i-1)*ny + j])*(old_C[(i+1)*ny + j] - old_C[(i-1)*ny + j]) +
										(old_v[i*ny + j] + old_v[i*ny + j-1])*(old_C[i*ny + j+1] - old_C[i*ny + j-1]));

			C[i*ny + j] +=  D*dt*hinv*hinv*(Ccopy[i*ny + j+1] - 4.*Ccopy[i*ny + j] + Ccopy[i*ny + j-1]
										  + Ccopy[(i+1)*ny + j] + Ccopy[(i-1)*ny + j]);
		}

	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
			old_C[i*ny + j] = Ccopy[i*ny + j];

	free(Ccopy);
}

void free_data_structure(Data *data)
{
	free(data->u);
	free(data->old_u);
	free(data->v);
	free(data->old_v);
	free(data->p);
	free(data->phi);
	free(data->T);
	free(data->old_T);
	free(data->Ca);
	free(data->old_Ca);
	free(data->Cb);
	free(data->old_Cb);
}
