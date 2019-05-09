/** * File              : ns_solver.c
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 25.04.2019
 * Last Modified Date: 09.05.2019
 */

#include <stdlib.h>
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
	double Tm = data->Tm;
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
			v[i*ny + j] -= 0.*beta*dt*(T[i*ny + j] - Tm)*g;
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
	unsigned long int ny = data->ny;

	/* double dt = data->dt; */
	/* double hinv = data->hinv; */

	double *u = data->u;
	double *v = data->v;
	/* double *old_v = data->old_v; */

	unsigned long int i;
	for (i = 1 + nx1; i <= 2*nx1; ++i) //Don't forget the useless ghost points.
		v[i*ny + ny-2] = u[i-3*nx1/4];
		/* v[i*ny + ny-2] = old_v[i*ny + ny-2] - hinv*dt*(old_v[i*ny + ny-2] - old_v[i*ny + ny-3]); */
		/* v[i*ny + ny-2] = old_v[i*ny + ny-2] - hinv*dt*(1.5*old_v[i*ny + ny-2] - 2*old_v[i*ny + ny-3] + .5*old_v[i*ny + ny-4]); */
		/* v[i*ny + ny-2] = v[i*ny + ny-3]; */
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

void free_data_structure(Data *data)
{
	free(data->u);
	free(data->old_u);
	free(data->v);
	free(data->old_v);
	free(data->p);
	free(data->phi);
	free(data->T);
	free(data->Ca);
	free(data->Cb);
}
