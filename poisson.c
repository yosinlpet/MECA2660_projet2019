/**
 * File              : poisson.c
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 06.05.2019
 * Last Modified Date: 09.05.2019
 */
/**File              : poisson.c
 * Date              : 26.04.2019
 * Last Modified Date: 06.05.2019
 */

#include <mpi.h>
#include "poisson.h"
#include "ns_solver.h"

/*Called by poisson_solver at each time step*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    -Impose zero mass flow here by changing value of U_star*/
void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, Data* data)
{
	unsigned long int nx = data->nx;
	unsigned long int ny = data->ny;
	double hinv = data->hinv;
	double dt = data->dt;
	double h = 1./hinv;

	double *u = data->u;
	double *v = data->v;

	double *divv = calloc(nx*ny, sizeof(double));

	unsigned long int i, j;
	for (i = 1; i < nx-1; ++i)
		for (j = 1; j < ny-1; ++j)
		{
			divv[i*ny + j] = u[i*ny + j] - u[(i-1)*ny + j] + v[i*ny + j] - v[i*ny + j-1];
		}
	
    int r;
	/* double a1 = 2*M_PI/3.; */
	/* double a2 = 2*M_PI/6.; */

	/* printf("===== RHS =====\n"); */
    for(r=rowStart+ny; r<rowEnd-ny ; r++)
	{
		if (is_boundary(r, nx, ny) == 0)
			rhs[r] = divv[r]*h/dt;
			/* rhs[r] = -(a1*a1*cos(a1*(r/ny - .5)*h) + a2*a2*cos(a2*(r % ny - .5)*h))*h*h; */
		else
			rhs[r] = 0.;
		/* printf("%f\t", rhs[r]); */
	}
	/* printf("\n===== END =====\n"); */

	free(divv);
}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
void poisson_solver(Poisson_data *pdata, Data* data)
{
	double *phi = data->phi;

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = pdata->sles;
    Vec b = pdata->b;
    Vec x = pdata->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(rhs, rowStart, rowEnd, data);
    VecRestoreArray(b, &rhs);

	/* PetscReal *uo = NULL; */
	/* VecNorm(b, NORM_2, uo); */
	/* printf("RHS NORM=%f\n", (double) *uo); */

    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    int r;
	/* printf("===== SOL =====\n"); */
    for(r=rowStart; r<rowEnd; r++)
	{
		phi[r] = sol[r];
		/* printf("%f\t", sol[r]); */
	}
	/* printf("\n===== END =====\n"); */

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
    /* free(uo); */
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
void computeLaplacianMatrix(Mat A, PetscInt rowStart, PetscInt rowEnd, unsigned long int nx, unsigned long int ny)
{
	int r;
	for (r = rowStart; r < rowEnd; r++)
	{
		/* MatSetValue(A, r, r, 1., INSERT_VALUES); */

		/* if (r == 38) */
		/* 	continue; */

		if(is_boundary(r, nx, ny) == 0)
		{
			MatSetValue(A, r, r, -4, INSERT_VALUES);
			MatSetValue(A, r, r+ny, 1, INSERT_VALUES);
			MatSetValue(A, r, r-ny, 1, INSERT_VALUES);
			MatSetValue(A, r, r+1, 1, INSERT_VALUES);
			MatSetValue(A, r, r-1, 1, INSERT_VALUES);
		}
		else if(is_boundary(r, nx, ny) == 1)
		{
			/* printf("RIGHT -- %d\n", r); */
			MatSetValue(A, r, r, 1, INSERT_VALUES);
			MatSetValue(A, r, r-ny, -1, INSERT_VALUES);
		}
		else if(is_boundary(r, nx, ny) == 2)
		{
			/* printf("LEFT -- %d\n", r); */
			MatSetValue(A, r, r, 1, INSERT_VALUES);
			MatSetValue(A, r, r+ny, -1, INSERT_VALUES);
		}
		else if(is_boundary(r, nx, ny) == 3)
		{
			/* printf("TOP -- %d\n", r); */
			MatSetValue(A, r, r, 1, INSERT_VALUES);
			MatSetValue(A, r, r-1, -1, INSERT_VALUES);
		}
		else
		{
			/* printf("BOT -- %d\n", r); */
			MatSetValue(A, r, r, 1, INSERT_VALUES);
			MatSetValue(A, r, r+1 , -1, INSERT_VALUES);
		}

	}
	/* MatSetValue(A, 38, 38, 1., INSERT_VALUES); */
}

int is_boundary(unsigned long int r, unsigned long int nx, unsigned long int ny)
{
	if(r >= ny*(nx-1))
		return 1; //right
	else if(r < ny)
		return 2; //left
	else if (r % ny == ny-1)
		return 3; //top
	else if(r % ny == 0)
		return 4; //bot
	else
		return 0;
}

/*To call during the initialization of your solver, before the begin of the time loop*/
PetscErrorCode initialize_poisson_solver(Poisson_data* data, unsigned long int nx, unsigned long int ny)
{
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	int nphi = nx*ny;

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A, 2*ny, NULL);
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(data->A, rowStart, rowEnd, nx, ny);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Singular matrix */
    MatNullSpace nullspace;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
    MatSetNullSpace(data->A, nullspace);
    MatSetTransposeNullSpace(data->A, nullspace);
    MatNullSpaceDestroy(&nullspace);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    PC prec;
    KSPGetPC(data->sles, &prec);
    KSPSetFromOptions(data->sles); // specifier le solveur exact a partir de la ligne commande mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-10, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,5);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Mattrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
void free_poisson_solver(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
