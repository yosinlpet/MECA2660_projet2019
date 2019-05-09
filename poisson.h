/**
 * File              : poisson.h
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 26.04.2019
 * Last Modified Date: 30.04.2019
 */
#ifndef _POISSON_H_
#define _POISSON_H_

/*To include in the file in which you will call initialize_poisson_solver and poisson_solver*/

#include <petsc.h>
#include <petscsys.h>
#include "ns_solver.h"

//Structure storing petsc vectors

typedef struct {

	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} Poisson_data;

PetscErrorCode initialize_poisson_solver(Poisson_data* data, unsigned long int nx, unsigned long int ny);
void poisson_solver(Poisson_data *pdata, Data *data);
void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, Data *data);
void computeLaplacianMatrix(Mat A, PetscInt rowStart, PetscInt rowEnd, unsigned long int nx, unsigned long int ny);
void free_poisson_solver(Poisson_data *data);
int is_boundary(unsigned long int r, unsigned long int nx, unsigned long int ny);

#endif

