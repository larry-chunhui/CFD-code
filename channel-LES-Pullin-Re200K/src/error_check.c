/*
 * serror_check.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "convective.h"
#include "share.h"
#include "derivatives__2.h"
#include "derivatives__4.h"
#include "interpolation.h"
#include "interpolation__2.h"
#include "interpolation__4.h"
#include "interpolation__les.h"
#include "error_check.h"
#include "les.h"
#include "spiral.h"
#include "velgrad_tensor.h"
#include "near_wall_ode.h"
#include "time_march.h"


void check_bc_dudy(variables *var0,
				   variables *var1,
				   int substep,
				   parameters *par)
{
	
	int i, k, index;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int nz = par->nz;
	const int local_nx_start = par->local_nx_start;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	
	const double gamma = par->rk3_gamma[substep];
	const double zeta = par->rk3_zeta[substep];
	const double dt_substep = (gamma+zeta)*par->dt;
	
	double *eta_true_t, *eta_true_b;
	eta_true_t = dvector(0, (local_nx+1)*(nz)-1);
	eta_true_b = dvector(0, (local_nx+1)*(nz)-1);
	
	
	double scratch;
	
	for (i=0; i<=(local_nx+1)*(nz)-1; ++i){
		scratch = 1.0/var0->eta0_b[i]*          exp( - dt_substep*var0->lambda_b[i]*var0->eta0_bar_b[i] )
        + 1.0/var0->eta0_bar_b[i]*(1.0- exp( - dt_substep*var0->lambda_b[i]*var0->eta0_bar_b[i] ));
		eta_true_b[i] = 1.0/scratch;
	}
	
	for (i=0; i<=(local_nx+1)*(nz)-1; ++i){
		scratch = 1.0/var0->eta0_t[i]*          exp( - dt_substep*var0->lambda_t[i]*var0->eta0_bar_t[i] )
        + 1.0/var0->eta0_bar_t[i]*(1.0- exp( - dt_substep*var0->lambda_t[i]*var0->eta0_bar_t[i] ));
		eta_true_t[i] = 1.0/scratch;
	}
    
	
	double error_t, error_b;
	double scale_t, scale_b;
	double local_error_t = 0.0;
	double local_error_b = 0.0;
	double local_scale_t = 0.0;
	double local_scale_b = 0.0;
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			index = par->index_tb[i][k];
			local_error_t += pow(eta_true_t[index] - var1->bc_dudy_t[index], 2.0);
			local_error_b += pow(eta_true_b[index] - var1->bc_dudy_b[index], 2.0);
			local_scale_t += pow(eta_true_t[index], 2.0);
			local_scale_b += pow(eta_true_b[index], 2.0);
		}
	}
	
	
	MPI_Reduce( &local_error_t, &error_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &local_error_b, &error_b, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &local_scale_t, &scale_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &local_scale_b, &scale_b, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (my_rank == 0){
		if (scale_t < 10e-16 || scale_b < 10e-16){
			error_t = sqrt(error_t)/(double)(par->nx*nz);
			error_b = sqrt(error_b)/(double)(par->nx*nz);
		}else{
			error_t = sqrt(error_t/scale_t);
			error_b = sqrt(error_b/scale_b);
		}
		printf("at substep %d wall ODE total error %e %e\n", substep, error_t, error_b);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
    
	free_dvector(eta_true_t, 0, (local_nx+1)*(nz)-1);
	free_dvector(eta_true_b, 0, (local_nx+1)*(nz)-1);
	
	
}



void check_wall_ode_time_march(
                               variables *var9,
                               variables *var0,
                               parameters *par,
                               int substep,
                               fftwplans *planptr)
{
	
	neighbors sh, *shared;
	shared = &sh;
	interpolated interp, *inter;
	inter = &interp;
	
	get_neignboring_data(var0, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
	get_periodic_data(var0, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
	init_interpolation(inter, par);
	interpolation__4(1.0, var0, shared, inter, par, planptr);
    
	MPI_Barrier(MPI_COMM_WORLD);
	
	get_convective(var0, shared, inter, par, planptr);
	
	get_dTijdxj__4(var0, inter, par, planptr);
	
	bc__nonlinear(var0, inter, substep, par, planptr);
	bc__rhs(var9, var0, substep, par);
    
	
	finalize_neignboring(shared, par);
	finalize_interpolation(inter, par);
	
}

