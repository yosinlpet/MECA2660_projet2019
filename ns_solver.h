/**
 * File              : ns_solver.h
 * Author            : Denis N Dumoulin <denis.dumoulin@uclouvain.be>
 * Date              : 26.04.2019
 * Last Modified Date: 07.05.2019
 */
#include "poisson.h"

#ifndef DEF_NS_SOLVER
#define DEF_NS_SOLVER

typedef struct {
	//parameters
	unsigned long int nx;
	unsigned long int ny;
	unsigned long int nx1;
	double hinv;
	double dt;
	double nu;
	double g;
	double Tm;
	double beta;
	
	//quantities that matter
	double *u;
	double *old_u;
	double *v;
	double *old_v;
	double *p;
	double *phi;
	
	//scalars
	double *T;
	double *Ca;
	double *Cb;
} Data;

void step_one(Data* data);
void reprojection(Data* data);
void impose_boundary(Data* data);
void flow_conservation(Data *data);
void free_data_structure(Data *data);
void ghost_points(Data *data, double vgl[]);
#endif
