/**
 * File              : projet.c
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 30.04.2019
 * Last Modified Date: 10.05.2019
 */
/**
 * File              : projet.c
 * Date              : 26.04.2019
 * Last Modified Date: 30.04.2019
 */
#include <stdlib.h> 
#include <math.h>
#include <petsc.h>
#include <petscmat.h>
#include <petscsys.h>

#include "poisson.h"
#include "ns_solver.h"
#include "Write.h"

static void show_usage() {
	printf("Usage: ./convection_diffusion <option(s)>\n\n");
	printf("Options:\n");
	printf("\t-h, --help\t\tShow this help message\n");
	printf("\t--nx\t\tNX\tNumber of cells on the width\n");
	printf("\t--cfl\t\tCFL\tCFL\n");
}

int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);

	double Re = 250.;
	double Pr = .7;
	double Sc = .1;
	double Gr = 2e5;
	double CFL = .05;
	unsigned long int nx1 = 64;
	double tinf = 50;

	for (int i=1; i < argc; ++i)
	{
		char *arg = argv[i];
		if (strcmp("-h", arg) == 0 || strcmp("--help", arg) == 0)
		{
			printf("[INFO] -- Help asked.\n");
			show_usage();
			return 1;
		}
		else if (strcmp("--nx", arg) == 0)
			if (i + 1 < argc)
				nx1 = atoi(argv[++i]);
			else
			{
				printf("[ERROR] -- NX requires an argument.\n");
				return 1;
			}
		else if (strcmp("--cfl", arg) == 0)
			if (i + 1 < argc)
				CFL = atof(argv[++i]);
			else
			{
				printf("[ERROR] -- CFL requires an argument.\n");
				return 1;
			}
		else
			break;
	}
	
	double hinv = (double) nx1;
	unsigned long int nx = 3*nx1 + 2;
	unsigned long int ny = 6*nx1 + 2;
	double nu = 1./Re;
	double alpha = nu/Pr;
	double D = nu/Sc;
	double g = -9.81;
	double beta = Gr*nu*nu/g;
	double k0 = 1e10*D;
	double ER = 20;
	double DH = 2*alpha/D;

	double dt = CFL/hinv;

	printf("CFL=%f\n", CFL);
	printf("dt=%f\n", dt);

	Data *data = malloc(sizeof(Data));
	data->nx1 = nx1;
	data->nx = nx;
	data->ny = ny;
	data->dt = dt;
	data->nu = nu;
	data->alpha = alpha;
	data->beta = beta;
	data->DH = DH;
	data->D = D;
	data->k0 = k0;
	data->ER = ER;
	data->hinv = hinv;

	Poisson_data *poisson = malloc(sizeof(Poisson_data));

	data->u = calloc(nx*ny, sizeof(double));
	data->v = calloc(nx*ny, sizeof(double));
	data->p = calloc(nx*ny, sizeof(double));
	data->phi = calloc(nx*ny, sizeof(double));
	data->old_u = calloc(nx*ny, sizeof(double));
	data->old_v = calloc(nx*ny, sizeof(double));

	data->T = calloc(nx*ny, sizeof(double));
	data->Ca = calloc(nx*ny, sizeof(double));
	data->Cb = calloc(nx*ny, sizeof(double));

	data->old_T = calloc(nx*ny, sizeof(double));
	data->old_Ca = calloc(nx*ny, sizeof(double));
	data->old_Cb = calloc(nx*ny, sizeof(double));

	double *vgl = calloc(ny, sizeof(double));
	double *u = data->u;
	double *v = data->v;
	double *T = data->T;
	double *Ca = data->Ca;
	double *old_u = data->old_u;
	double *old_v = data->old_v;
	double *old_T = data->old_T;
	double *old_Ca = data->old_Ca;

	unsigned long int i, j;
	//Initialize atmosphere
	for (i = 0; i < nx; ++i)
		for (j = 0; j < ny; ++j)
		{
			T[i*ny + j] = .6;
			Ca[i*ny + j] = .21;
		}

	//Initialize flow at the inlet
	for (j = nx1/4 + 1; j <= 5*nx1/4; ++j)
	{
		vgl[j] = -1.5*(1. - 4*pow(((double) (j-1))/hinv - .75, 2.))*sin(.25*M_PI); 
		u[j] = 1.5*(1. - 4*pow(((double) (j-1) + .5)/hinv - .75, 2.))*cos(.25*M_PI); 
		T[j]  = 1.;
	}
	impose_boundary(data);	
	flow_conservation(data);
	ghost_points(data, vgl);

	for (i = 0; i < nx; ++i)
	{
		for(j = 0; j < nx; ++j)
		{
			old_u[i*ny + j] = u[i*ny + j];
			old_v[i*ny + j] = v[i*ny + j];
			old_T[i*ny + j] = T[i*ny + j];
			old_Ca[i*ny + j] = Ca[i*ny + j];
		}
	}

	initialize_poisson_solver(poisson, nx, ny);

	//Time loop
	unsigned long int n = 0;
	double t;
	for (t = 0.; t <= tinf; t += dt) 
	{
		printf("\n%f/%f\n", t, tinf);
		step_one(data);
		flow_conservation(data);
		poisson_solver(poisson, data);
		reprojection(data);
		if (n % 100 == 0)
			logger(data, n);

		temp(data);
		concentration(data, 'a');
		concentration(data, 'b');

		impose_boundary(data);	
		ghost_points(data, vgl);
		++n;
	}

	//Free memory
	free(vgl);
	free_poisson_solver(poisson);
	free(poisson);
    PetscFinalize();
	free_data_structure(data);
	free(data);
}
