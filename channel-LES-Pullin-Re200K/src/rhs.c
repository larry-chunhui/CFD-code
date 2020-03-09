/*
 *  rhs.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Calculates the right hand side of the first equation (of the three). 
 *          - void get_rhs_u/v/w(variables *var1, variables *var0, variables *var9, double *rhs_u, int substep, parameters *par)
 * Calculates the right hand side of the pressure poisson equation (second of the three).
 *          - void get_rhs_poisson(variables *var, double *rhs_poisson, int substep, parameters *par)
 * Please refer to the document associated with this code for the equations.
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "definitions.h"
#include "nrutil.h"
#include "share.h"
#include "derivatives__2.h"
#include "derivatives__4.h"
#include "interpolation.h"
#include "interpolation__2.h"
#include "interpolation__4.h"
#include "rhs.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

/* ******************************************************************* */

/* ******************************************************************* */
void get_rhs_u(
               variables *var1,
               variables *var0,
               variables *var9,
               double *rhs_u,
               int substep,
               parameters *par)
{
    
	int i;
	const int local_size_u = par->local_size_u;
	const double dt = par->dt;
	const double re = par->re;
	
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
	const double gamma = par->rk3_gamma[substep];
	const double zeta = par->rk3_zeta[substep];
	
	
	for (i=0; i<=local_size_u-1; ++i){
		rhs_u[i] = var0->u[i] - dt*(alpha+beta)*var0->Gpx[i] \
        + dt*alpha*var0->Lu_y[i]/re \
        - dt*(gamma*var0->Nu[i] + zeta*var9->Nu[i]) \
        + dt*(gamma*var0->Lu_xz[i] + zeta*var9->Lu_xz[i])/re \
        + dt*beta*var1->bc_Lu_y[i]/re;
	}
	
	return;
}

/* ******************************************************************* */

/* ******************************************************************* */
void get_rhs_v(
               variables *var1,
               variables *var0,
               variables *var9,
               double *rhs_v,
               int substep,
               parameters *par)
{
	int i;
	const int local_size_v = par->local_size_v;
	const double dt = par->dt;
	const double re = par->re;
	
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
	const double gamma = par->rk3_gamma[substep];
	const double zeta = par->rk3_zeta[substep];
	
	for (i=0; i<=local_size_v-1; ++i){
		rhs_v[i] = var0->v[i] - dt*(alpha+beta)*var0->Gpy[i]  \
        + dt*alpha*var0->Lv_y[i]/re \
        - dt*(gamma*var0->Nv[i] + zeta*var9->Nv[i]) \
        + dt*(gamma*var0->Lv_xz[i] + zeta*var9->Lv_xz[i])/re \
        + dt*beta*var1->bc_Lv_y[i]/re;
	}
    
	
	return;
}

/* ******************************************************************* */

/* ******************************************************************* */
void get_rhs_w(
               variables *var1,
               variables *var0,
               variables *var9,
               double *rhs_w,
               int substep,
               parameters *par)
{
	int i;
	const int local_size_wp = par->local_size_wp;
	const double dt = par->dt;
	const double re = par->re;
	
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
	const double gamma = par->rk3_gamma[substep];
	const double zeta = par->rk3_zeta[substep];
	
	for (i=0; i<=local_size_wp-1; ++i){
		rhs_w[i] = var0->w[i] - dt*(alpha+beta)*var0->Gpz[i] \
        + dt*alpha*var0->Lw_y[i]/re \
        - dt*(gamma*var0->Nw[i] + zeta*var9->Nw[i]) \
        + dt*(gamma*var0->Lw_xz[i] + zeta*var9->Lw_xz[i])/re \
        + dt*beta*var1->bc_Lw_y[i]/re;
	}
	
	
	
	
	return;
}


/* ******************************************************************* */

/* ******************************************************************* */
void get_rhs_poisson(variables *var,
					 double *rhs_poisson,
					 int substep,
					 parameters *par)
{
    
	neighbors sh, *shared;
	shared = &sh;
	interpolated interp, *inter;
	inter = &interp;
	
	int i;
	const int local_size_wp = par->local_size_wp;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
	
	get_neignboring_data(var, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
    
	get_periodic_data(var, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
	inter->u_ued = FFTW_MALLOC(local_size_ued);
	inter->v_ved = FFTW_MALLOC(local_size_ved);
	
	
	if (par->fd_order == 4){
		bigger_array_u_ued__4(1.0, var, shared, inter->u_ued, par);
		bigger_array_v_ved__4(1.0, var, shared, inter->v_ved, par);
		MPI_Barrier(MPI_COMM_WORLD);
		get_divergence__4(inter->u_ued, inter->v_ved, var->w, var->Du, par);
	}
	if (par->fd_order == 2){
		bigger_array_u_ued__2(1.0, var, shared, inter->u_ued, par);
		bigger_array_v_ved__2(1.0, var, shared, inter->v_ved, par);
		MPI_Barrier(MPI_COMM_WORLD);
		get_divergence__2(inter->u_ued, inter->v_ved, var->w, var->Du, par);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	finalize_neignboring(shared, par);
	fftw_free(inter->u_ued);
	fftw_free(inter->v_ved);
	
	for (i=0; i<=local_size_wp-1; ++i){
		rhs_poisson[i] = var->Du[i]/(alpha+beta)/par->dt;
	}
	
	
	return;
}

