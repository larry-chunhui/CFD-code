/*
 *  time_march.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * time_march() runs at each fractional step with different values for
 * alpha, beta, gamma and zeta (coefficients for the numerical method). 
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "convective.h"
#include "share.h"
#include "statistics.h"
#include "differentiate.h"
#include "derivatives__2.h"
#include "derivatives__4.h"
#include "interpolation.h"
#include "interpolation__2.h"
#include "interpolation__4.h"
#include "interpolation__les.h"
#include "trisolver.h"
#include "profiles.h"
#include "pentasolver.h"
#include "poisson_solver_transpose.h"
#include "correction.h"
#include "rhs.h"
#include "velgrad_tensor.h"
#include "bc__get.h"
#include "near_wall_ode.h"
#include "les.h"
#include "error_check.h"
#include "time_march.h"

extern void fill_tb(double *bc, parameters *par, int tag);

/* ******************************************************************* */

/* ******************************************************************* */
void time_march(
				variables *var9,
				variables *var0,
				variables *var1,
				parameters *par,
				int substep,
				fftwplans *planptr)
{
	
	neighbors sh, *shared;
	shared = &sh;
	interpolated interp, *inter;
	inter = &interp;
	
	int i;
	const int local_nx = par->local_nx;
	const int nz = par->nz;
	const int local_size_u = par->local_size_u;
	const int local_size_v = par->local_size_v;
	const int local_size_wp = par->local_size_wp;
    
	double *rhs_u;
	double *rhs_v;
	double *rhs_w;
	double *rhs_poisson;
	
	
	
	rhs_u = dvector(0, local_size_u-1); // freed
	rhs_v = dvector(0, local_size_v-1); // freed
	rhs_w = dvector(0, local_size_wp-1); // freed
	rhs_poisson = dvector(0, local_size_wp-1); // freed
	
	
	
	get_neignboring_data(var0, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
    
	get_periodic_data(var0, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
	init_interpolation(inter, par);
	
//if(par->my_rank == 0){
//printf("time march/ init_interpolation succeed\n");
//}	
	
	if (par->fd_order == 2){
		interpolation__2(1.0, var0, shared, inter, par, planptr);
	}
	if (par->fd_order == 4){
		interpolation__4(1.0, var0, shared, inter, par, planptr);
	}
	MPI_Barrier(MPI_COMM_WORLD);

//if(par->my_rank == 0){
//printf("time march/ interpolation__4 succeed\n");
//}		
	
	
	get_convective(var0, shared, inter, par, planptr);
//if(par->my_rank == 0){
//printf("time march/ get_convective succeed\n");
//}	
	
	if (par->fd_order == 2){
		get_dTijdxj__2(var0, inter, par, planptr);
	}
	if (par->fd_order == 4){
		get_dTijdxj__4(var0, inter, par, planptr);
	}
//if(par->my_rank == 0){
//printf("channel integration/time march/ get_dTijdxj__4 succeed\n");      //已成功
//}	
	
	if (par->bc == 1){
		bc__nonlinear(var0, inter, substep, par, planptr);
		bc__rhs(var9, var0, substep, par);
		bc__time_march_bc_dudy(var9, var0, var1, par);
		
		term_check__nonlinear_ode(var0, inter, substep, par, planptr);
		//bc__time_march_bc_dudy_analitic(var0, var1, substep, par);
		check_bc_dudy(var0, var1, substep, par);
		
		fftw_execute_r2r(planptr->p1d_z_bc_tb, var1->bc_dudy_t, var1->bc_dudy_t);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, var1->bc_dudy_b, var1->bc_dudy_b);
		for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
			var1->bc_dudy_t[i] = var1->bc_dudy_t[i]/(double)(nz);
			var1->bc_dudy_b[i] = var1->bc_dudy_b[i]/(double)(nz);
		}
		bc_trancate_exp(var1->bc_dudy_t, "tb", par);
		bc_trancate_exp(var1->bc_dudy_b, "tb", par);
		
		fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_dudy_t, var1->bc_dudy_t);
		fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_dudy_b, var1->bc_dudy_b);
		
		
		bc__get_bc_vel(var1->bc_dudy_t, var0->bc_kappa_t, var0->bc_kappa_t,
					   var1->bc_u_t, var1->bc_v_t, var1->bc_w_t, "top", par);
		
		bc__get_bc_vel(var1->bc_dudy_b, var0->bc_kappa_b, var0->bc_kappa_b,
					   var1->bc_u_b, var1->bc_v_b, var1->bc_w_b, "bot", par);
        
		
		fftw_execute_r2r(planptr->p1d_z_bc_tb, var1->bc_u_t, var1->bc_u_t);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, var1->bc_u_b, var1->bc_u_b);
		for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
			var1->bc_u_t[i] = var1->bc_u_t[i]/(double)(nz);
			var1->bc_u_b[i] = var1->bc_u_b[i]/(double)(nz);
		}
		
		bc_trancate_exp(var1->bc_u_t, "tb", par);
		bc_trancate_exp(var1->bc_u_b, "tb", par);
		
        
        const int tag = 223;
        fill_tb(var1->bc_u_t, par, tag);
        fill_tb(var1->bc_v_t, par, tag+1);
        fill_tb(var1->bc_w_t, par, tag+2);
        
        fill_tb(var1->bc_u_b, par, tag+3);
        fill_tb(var1->bc_v_b, par, tag+4);
        fill_tb(var1->bc_w_b, par, tag+5);
        
	}
	
//if(par->my_rank == 0){
//printf("channel integration/time march/ fill_tb succeed\n"); //已成功
//}    
	
	
	if (par->fd_order == 4){
		get_laplace_bc_y__4(var1, par);
		get_laplace_y__4(inter->u_ued, var0->Lu_y, inter->v_ved,
                         var0->Lv_y, inter->w_ct, var0->Lw_y, var0, par);
		get_laplace_xz__4(inter->u_ued, var0->Lu_xz, inter->v_ved,
						  var0->Lv_xz, inter->w_ct, var0->Lw_xz, var0, par);
	}
	if (par->fd_order == 2){
		get_laplace_bc_y__2(var1, par);
		get_laplace_y__2(inter->u_ued, var0->Lu_y, inter->v_ved,
						 var0->Lv_y, inter->w_ct, var0->Lw_y, var0, par);
		get_laplace_xz__2(inter->u_ued, var0->Lu_xz, inter->v_ved,
						  var0->Lv_xz, inter->w_ct, var0->Lw_xz, var0, par);
	}
	MPI_Barrier(MPI_COMM_WORLD);

//if(par->my_rank == 0){
//printf("channel integration/time march/get_laplace succeed\n"); //已成功
//}  	
	
	
//if(par->my_rank == 0){
//printf("channel integration/time march/ statistics _les begin\n"); 
//printf("par->it=%d,par->stat=%d,substep=%d,par->les=%d\n",par->it,par->stat,substep,par->les); 
//}  	
	if (par->it%par->stat == 0 && substep == 0){
	 
		statistics(inter, planptr, par);       //此地有问题
	
	}

	if (par->it%par->stat == 0 && substep == 0 && par->les == 1){
	
		statistics_les(inter, planptr, par);
			
	}


	finalize_neignboring(shared, par);
	finalize_interpolation(inter, par);
	
	
    
	get_rhs(var1, var0, var9, rhs_u, rhs_v, rhs_w, substep, par);
	MPI_Barrier(MPI_COMM_WORLD);

	
	
	
	
	
	// solve for u star
	vel_solvers(var1->u, var1->v, var1->w, rhs_u, rhs_v, rhs_w, substep, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
    
	
	get_rhs_poisson(var1, rhs_poisson, substep, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	poisson_1DDCT_in_y_transpose_1DFFT_in_x(var1->p, rhs_poisson, par, planptr);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	get_correction(var1, substep, par);
	MPI_Barrier(MPI_COMM_WORLD);
		
	
	// update field variables
	for (i=0; i<=local_size_u-1; ++i){
		var9->Nu[i] = var0->Nu[i];
		var9->Lu_xz[i] = var0->Lu_xz[i];
		var0->bc_Lu_y[i] = var1->bc_Lu_y[i];
		var0->u[i] = var1->u[i];
		
		var0->Gpx[i] = var0->Gpx[i] + var1->Gpx[i];
	}
	for (i=0; i<=local_size_v-1; ++i){
		var9->Nv[i] = var0->Nv[i];
		var9->Lv_xz[i] = var0->Lv_xz[i];
		var0->bc_Lv_y[i] = var1->bc_Lv_y[i];
		var0->v[i] = var1->v[i];
		
		var0->Gpy[i] = var0->Gpy[i] + var1->Gpy[i];
	}
	for (i=0; i<=local_size_wp-1; ++i){
		var9->Nw[i] = var0->Nw[i];
		var9->Lw_xz[i] = var0->Lw_xz[i];
		var0->bc_Lw_y[i] = var1->bc_Lw_y[i];
		var0->w[i] = var1->w[i];
		
		var0->Gpz[i] = var0->Gpz[i] + var1->Gpz[i];
		var0->p[i] = var0->p[i] + var1->p[i];
	}
	
	
	
	if (par->les == 1 && par->bc == 1){
		for (i=0; i<=par->local_size_tb-1; ++i){
			var0->bc_dudy_t[i] = var1->bc_dudy_t[i];
			var0->bc_u_t[i] = var1->bc_u_t[i];
			var0->bc_v_t[i] = var1->bc_v_t[i];
			var0->bc_w_t[i] = var1->bc_w_t[i];
			
			
			var0->bc_dudy_b[i] = var1->bc_dudy_b[i];
			var0->bc_u_b[i] = var1->bc_u_b[i];
			var0->bc_v_b[i] = var1->bc_v_b[i];
			var0->bc_w_b[i] = var1->bc_w_b[i];
			
		}
	}
if(par->my_rank == 0){
// printf("update velocity field succeed\n");
} 	
	free_dvector(rhs_u, 0, local_size_u-1);
	free_dvector(rhs_v, 0, local_size_v-1);
	free_dvector(rhs_w, 0, local_size_wp-1);
	free_dvector(rhs_poisson, 0, local_size_wp-1);
	
	
	return;
}


/* ******************************************************************* */

/* ******************************************************************* */
void get_rhs(
             variables *var1,
             variables *var0,
             variables *var9,
             double *rhs_u,
             double *rhs_v,
             double *rhs_w,
             int substep,
             parameters *par)
{
    
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
    
	get_rhs_u(var1, var0, var9, rhs_u, substep, par);
	get_rhs_v(var1, var0, var9, rhs_v, substep, par);
	get_rhs_w(var1, var0, var9, rhs_w, substep, par);
    
	// cosntant pressure gradient
	if (par->const_mass == 0){
		int i, j, k;
		const int start = par->local_nx_start;
		const int end = par->local_nx_start + par->local_nx;
		const double dt = par->dt;
		
		k = 0;
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=par->ny-1; ++j){
				rhs_u[par->index_u[i][j][k]] += - dt*(alpha+beta)*par->dPdx;
			}
		}
	}
	
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void get_rhs_transpose(
                       variables *var1,
                       variables *var0,
                       variables *var9,
                       double *rhs_u,
                       double *rhs_v,
                       double *rhs_w,
                       int substep,
                       parameters *par)
{
    
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
    
	get_rhs_u(var1, var0, var9, rhs_u, substep, par);
	get_rhs_v(var1, var0, var9, rhs_v, substep, par);
	get_rhs_w(var1, var0, var9, rhs_w, substep, par);
    
	// cosntant pressure gradient
	if (par->const_mass == 0){
		int i, j, k;
		const int start = par->local_nx_start;
		const int end = par->local_nx_start + par->local_nx;
		const double dt = par->dt;
		
		k = 0;
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=par->ny-1; ++j){
				rhs_w[par->index_wp[i][j][k]] += - dt*(alpha+beta)*par->dPdx;
			}
		}
		
		
		
	}
	
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void vel_solvers(
                 double *u,
                 double *v,
                 double *w,
                 double *rhs_u,
                 double *rhs_v,
                 double *rhs_w,
                 int substep,
                 parameters *par)
{
	
	int i, j;
	const int my_rank = par->my_rank;
	const int ny = par->ny;
	const int local_size_u = par->local_size_u;
	const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
    
	
	
    
	if (par->fd_order == 2){
		vel_trisolver_v(v, rhs_v, substep, par);
		vel_trisolver_w(w, rhs_w, substep, par);
	}
	
	if (par->fd_order == 4){
		vel_pentasolver_v(v, rhs_v, substep, par);
		vel_pentasolver_w(w, rhs_w, substep, par);
	}
    
	
    
    
	if (par->const_mass == 0){
		// const pressure gradient case
		if (par->fd_order == 2){
			vel_trisolver_u(u, rhs_u, substep, par);
		}
		
		if (par->fd_order == 4){
			vel_pentasolver_u(u, rhs_u, substep, par);
		}
        
        
	}else{
		// constant mass flux case
		double *u1, *u2;
		double *rhs;
		double coeff;
		double local_bulk_u1, local_bulk_u2;
		double bulk_u1, bulk_u2;
		
		rhs = dvector(0, local_size_u-1); // freed
		u1 = dvector(0, local_size_u-1); // freed
		u2 = dvector(0, local_size_u-1); // freed
		
		for (i=0; i<=local_size_u-1; ++i){
			rhs[i] = 0.0;
			u1[i] = 0.0;
			u2[i] = 0.0;
		}
		
		
		if (par->fd_order == 2){
			vel_trisolver_u(u1, rhs_u, substep, par);
		}
		
		if (par->fd_order == 4){
			vel_pentasolver_u(u1, rhs_u, substep, par);
		}
		
		
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				rhs[par->index_u[i][j][0]] = - 1.0;
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (par->fd_order == 2){
			vel_trisolver_u(u2, rhs, substep, par);
		}
		if (par->fd_order == 4){
			vel_pentasolver_u(u2, rhs, substep, par);
		}
		
		MPI_Bcast(u2, local_size_u, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
        // compute mass flux (bulk_u) with u = u1 + u2
		bulk_u1 = 0.0;
		bulk_u2 = 0.0;
		local_bulk_u1 = 0.0;
		local_bulk_u2 = 0.0;
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				local_bulk_u1 += u1[par->index_u[i][j][0]]/(double)(local_nx*ny);
				local_bulk_u2 += u2[par->index_u[i][j][0]]/(double)(local_nx*ny);
			}
		}
		
		MPI_Reduce(&local_bulk_u1, &bulk_u1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&local_bulk_u2, &bulk_u2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (my_rank == 0){
			bulk_u1 = bulk_u1/(double)par->nproc;
			bulk_u2 = bulk_u2/(double)par->nproc;
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Bcast(&bulk_u1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&bulk_u2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		// compute alpha for u = u1 + alpha*u2 to satisfy const par->mass_flux
		coeff = (par->mass_flux - bulk_u1)/bulk_u2;
		par->dPdx = coeff/par->dt/(alpha+beta);
		
		// compute u = u1 + alpha*u2
		for (i=0; i<=local_size_u-1; ++i){
			u[i] = u1[i] + coeff*u2[i];
			rhs_u[i] += - 1.0*coeff;
		}
		
		free_dvector(rhs, 0, local_size_u-1);
		free_dvector(u1, 0, local_size_u-1);
		free_dvector(u2, 0, local_size_u-1);
		
    }
    
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void get_dTijdxj__4(
                    variables *var,
                    interpolated *inter,
                    parameters *par,
                    fftwplans *planptr)
{
	
	int i;
	double *scratch;
	
	if (par->les == 0) {
		// set all subgrid stress to be zero
		set_Tij_zero(var, inter, par);
	}
	
	if (par->les == 1){
		// get velocity gradient (physical space)
		get_velgrad_tensor_cn(var, inter, par, planptr);
		// get subgrid stress (Fourier space)
		les_channel(var, inter, par, planptr);
		// interpolate subgrid stress (Fourier space)
		interpolation_les(var, inter, par);
        
        
		// ************************************************
		// ************************************************
        
        
		scratch = FFTW_MALLOC(par->local_size_u); // freed
		// x-components
		// dTxx/dx, dTxy/dy, dTzx/dz
        
		// ************************************************
		// dTxx/dx
		derivative_x1_ct2u(inter->Txx_ct, scratch, par);
        
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += 9.0/8.0*scratch[i];
		}
		derivative_x3_ct2u(inter->Txx_ct, scratch, par);
        
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += - 1.0/8.0*scratch[i];
		}
        
		// ************************************************
		// dTxy/dy
		derivative_y1_cn2u(inter->Txy_cn, scratch, par);
		
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += 9.0/8.0*scratch[i];
		}
        
		derivative_y3_cn2u(inter->Txy_cn, scratch, par);
		
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += - 1.0/8.0*scratch[i];
		}
        
		// ************************************************
		// dTzx/dz
		derivative_z_ued2u(inter->Tzx_ued, scratch, par);
		
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += scratch[i];
		}
        
		fftw_free(scratch);
        
        
        
        
		
		// ************************************************
		// ************************************************
		scratch = FFTW_MALLOC(par->local_size_v); // freed
		// y-components
		// dTxy/dx, dTyy/dy, dTyz/dz
        
		// ************************************************
		// dTxy/dx
		derivative_x1_cn2v(inter->Txy_cn, scratch, par);
        
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += 9.0/8.0*scratch[i];
		}
		derivative_x3_cn2v(inter->Txy_cn, scratch, par);
        
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += - 1.0/8.0*scratch[i];
		}
        
		// ************************************************
		// dTyy/dy
		derivative_y1_ct2v(inter->Tyy_ct, scratch, par);
        
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += 9.0/8.0*scratch[i];
		}
        
		derivative_y3_ct2v(inter->Tyy_ct, scratch, par);
        
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += - 1.0/8.0*scratch[i];
		}
        
		// ************************************************
		// dTyz/dz
		derivative_z_ved2v(inter->Tyz_ved, scratch, par);
        
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += scratch[i];
		}
        
		fftw_free(scratch);
        
        
        
        
		// ************************************************
		// ************************************************
		scratch = FFTW_MALLOC(par->local_size_wp); // freed
		// z-components
		// dTzx/dx, dTyz/dy, dTzz/dz
        
		// ************************************************
		// dTzx/dx
		derivative_x1_ued2wp(inter->Tzx_ued, scratch, par);
        
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += 9.0/8.0*scratch[i];
		}
		derivative_x3_ued2wp(inter->Tzx_ued, scratch, par);
        
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += - 1.0/8.0*scratch[i];
		}
        
		// ************************************************
		// dTyz/dy
		derivative_y1_ved2wp(inter->Tyz_ved, scratch, par);
        
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += 9.0/8.0*scratch[i];
		}
        
		derivative_y3_ved2wp(inter->Tyz_ved, scratch, par);
        
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += - 1.0/8.0*scratch[i];
		}
        
		// ************************************************
		// dTzz/dz
		derivative_z_ct2wp(inter->Tzz_ct, scratch, par);
        
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += scratch[i];
		}
        
		fftw_free(scratch);
        
	}
	
	return;
}


/* ******************************************************************* */

/* ******************************************************************* */
void get_dTijdxj__2(
					variables *var,
					interpolated *inter,
					parameters *par,
					fftwplans *planptr)
{
	
	int i;
	double *scratch;
	
	if (par->les == 0) {
		// set all subgrid stress to be zero
		set_Tij_zero(var, inter, par);
	}
	
	if (par->les == 1){
        
        
		// get velocity gradient (physical space)
		get_velgrad_tensor_cn__2(var, inter, par, planptr);
		// get subgrid stress (Fourier space)
		les_channel(var, inter, par, planptr);
		// interpolate subgrid stress (Fourier space)
		interpolation_les__2(var, inter, par);
		
		
		// ************************************************
		// ************************************************
		
		
		scratch = FFTW_MALLOC(par->local_size_u); // freed
		// x-components
		// dTxx/dx, dTxy/dy, dTzx/dz
		
		// ************************************************
		// dTxx/dx
		derivative_x1_ct2u(inter->Txx_ct, scratch, par);
		
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += scratch[i];
		}
		
		// ************************************************
		// dTxy/dy
		derivative_y1_cn2u(inter->Txy_cn, scratch, par);
		
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += scratch[i];
		}
		
		
		// ************************************************
		// dTzx/dz
		derivative_z_ued2u(inter->Tzx_ued, scratch, par);
		
		for (i=0; i<=par->local_size_u-1; ++i){
			var->Nu[i] += scratch[i];
		}
		
		fftw_free(scratch);
		
		
		
		
		
		// ************************************************
		// ************************************************
		scratch = FFTW_MALLOC(par->local_size_v); // freed
		// y-components
		// dTxy/dx, dTyy/dy, dTyz/dz
		
		// ************************************************
		// dTxy/dx
		derivative_x1_cn2v(inter->Txy_cn, scratch, par);
		
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += scratch[i];
		}
		
		// ************************************************
		// dTyy/dy
		derivative_y1_ct2v(inter->Tyy_ct, scratch, par);
		
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += scratch[i];
		}
		
		
		// ************************************************
		// dTyz/dz
		derivative_z_ved2v(inter->Tyz_ved, scratch, par);
		
		for (i=0; i<=par->local_size_v-1; ++i){
			var->Nv[i] += scratch[i];
		}
		
		fftw_free(scratch);
		
		
		
		
		// ************************************************
		// ************************************************
		scratch = FFTW_MALLOC(par->local_size_wp); // freed
		// z-components
		// dTzx/dx, dTyz/dy, dTzz/dz
		
		// ************************************************
		// dTzx/dx
		derivative_x1_ued2wp(inter->Tzx_ued, scratch, par);
		
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += scratch[i];
		}
		
		// ************************************************
		// dTyz/dy
		derivative_y1_ved2wp(inter->Tyz_ved, scratch, par);
		
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += scratch[i];
		}
		
		
		// ************************************************
		// dTzz/dz
		derivative_z_ct2wp(inter->Tzz_ct, scratch, par);
		
		for (i=0; i<=par->local_size_wp-1; ++i){
			var->Nw[i] += scratch[i];
		}
		
		fftw_free(scratch);
		
	}
	
	return;
}


/* ******************************************************************* */
// 4th order test
/* ******************************************************************* */
void get_dTijdxj__test(
                       variables *var,
                       interpolated *inter,
                       double *dTxxdx,
                       double *dTxydy,
                       double *dTzxdz,
                       double *dTxydx,
                       double *dTyydy,
                       double *dTyzdz,
                       double *dTzxdx,
                       double *dTyzdy,
                       double *dTzzdz,
                       parameters *par,
                       fftwplans *planptr)
{
	
	int i;
	double *scratch;
	
	
	// ************************************************
	// ************************************************
	scratch = FFTW_MALLOC(par->local_size_u);
	// x-components
	// dTxx/dx, dTxy/dy, dTzx/dz
	
	// ************************************************
	// dTxx/dx
	derivative_x1_ct2u(inter->Txx_ct, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTxxdx[i] = 9.0/8.0*scratch[i];
	}
	derivative_x3_ct2u(inter->Txx_ct, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTxxdx[i] += - 1.0/8.0*scratch[i];
	}
	
	// ************************************************
	// dTxy/dy
	derivative_y1_cn2u(inter->Txy_cn, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTxydy[i] = 9.0/8.0*scratch[i];
	}
	
	derivative_y3_cn2u(inter->Txy_cn, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTxydy[i] += - 1.0/8.0*scratch[i];
	}
	
	// ************************************************
	// dTzx/dz
	derivative_z_ued2u(inter->Tzx_ued, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTzxdz[i] = scratch[i];
	}
	
	fftw_free(scratch);
	
	
	
	// ************************************************
	// ************************************************
	scratch = FFTW_MALLOC(par->local_size_v);
	// y-components
	// dTxy/dx, dTyy/dy, dTyz/dz
	
	// ************************************************
	// dTxy/dx
	derivative_x1_cn2v(inter->Txy_cn, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTxydx[i] = 9.0/8.0*scratch[i];
	}
	derivative_x3_cn2v(inter->Txy_cn, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTxydx[i] += - 1.0/8.0*scratch[i];
	}
	
	// ************************************************
	// dTyy/dy
	derivative_y1_ct2v(inter->Tyy_ct, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTyydy[i] = 9.0/8.0*scratch[i];
	}
	
	derivative_y3_ct2v(inter->Tyy_ct, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTyydy[i] += - 1.0/8.0*scratch[i];
	}
	
	// ************************************************
	// dTyz/dz
	derivative_z_ved2v(inter->Tyz_ved, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTyzdz[i] = scratch[i];
	}
	
	fftw_free(scratch);
	
	
	// ************************************************
	// ************************************************
	scratch = FFTW_MALLOC(par->local_size_wp);
	// z-components
	// dTzx/dx, dTyz/dy, dTzz/dz
	
	// ************************************************
	// dTzx/dx
	derivative_x1_ued2wp(inter->Tzx_ued, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTzxdx[i] = 9.0/8.0*scratch[i];
	}
	derivative_x3_ued2wp(inter->Tzx_ued, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTzxdx[i] += - 1.0/8.0*scratch[i];
	}
	
	// ************************************************
	// dTyz/dy
	derivative_y1_ved2wp(inter->Tyz_ved, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTyzdy[i] = 9.0/8.0*scratch[i];
	}
	
	derivative_y3_ved2wp(inter->Tyz_ved, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTyzdy[i] += - 1.0/8.0*scratch[i];
	}
	
	// ************************************************
	// dTzz/dz
	derivative_z_ct2wp(inter->Tzz_ct, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTzzdz[i] = scratch[i];
	}
	
	fftw_free(scratch);
	
	
	return;
}



/* ******************************************************************* */
// 2nd order test
/* ******************************************************************* */
void get_dTijdxj__2_test(
                         variables *var,
                         interpolated *inter,
                         double *dTxxdx,
                         double *dTxydy,
                         double *dTzxdz,
                         double *dTxydx,
                         double *dTyydy,
                         double *dTyzdz,
                         double *dTzxdx,
                         double *dTyzdy,
                         double *dTzzdz,
                         parameters *par,
                         fftwplans *planptr)
{
	
	int i;
	double *scratch;
    
    
	scratch = FFTW_MALLOC(par->local_size_u); // freed
	// x-components
	// dTxx/dx, dTxy/dy, dTzx/dz
	
	// ************************************************
	// dTxx/dx
	derivative_x1_ct2u(inter->Txx_ct, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTxxdx[i] = scratch[i];
	}
	
	// ************************************************
	// dTxy/dy
	derivative_y1_cn2u(inter->Txy_cn, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTxydy[i] = scratch[i];
	}
	
	
	// ************************************************
	// dTzx/dz
	derivative_z_ued2u(inter->Tzx_ued, scratch, par);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		dTzxdz[i] = scratch[i];
	}
	
	fftw_free(scratch);
	
	
	
	
	
	// ************************************************
	// ************************************************
	scratch = FFTW_MALLOC(par->local_size_v); // freed
	// y-components
	// dTxy/dx, dTyy/dy, dTyz/dz
	
	// ************************************************
	// dTxy/dx
	derivative_x1_cn2v(inter->Txy_cn, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTxydx[i] = scratch[i];
	}
	
	// ************************************************
	// dTyy/dy
	derivative_y1_ct2v(inter->Tyy_ct, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTyydy[i] = scratch[i];
	}
	
	
	// ************************************************
	// dTyz/dz
	derivative_z_ved2v(inter->Tyz_ved, scratch, par);
	
	for (i=0; i<=par->local_size_v-1; ++i){
		dTyzdz[i] = scratch[i];
	}
	
	fftw_free(scratch);
	
	
	
	
	// ************************************************
	// ************************************************
	scratch = FFTW_MALLOC(par->local_size_wp); // freed
	// z-components
	// dTzx/dx, dTyz/dy, dTzz/dz
	
	// ************************************************
	// dTzx/dx
	derivative_x1_ued2wp(inter->Tzx_ued, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTzxdx[i] = scratch[i];
	}
	
	// ************************************************
	// dTyz/dy
	derivative_y1_ved2wp(inter->Tyz_ved, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTyzdy[i] = scratch[i];
	}
	
	
	// ************************************************
	// dTzz/dz
	derivative_z_ct2wp(inter->Tzz_ct, scratch, par);
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		dTzzdz[i] = scratch[i];
	}
	
	fftw_free(scratch);
	
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void set_Tij_zero(variables *var, interpolated *inter, parameters *par)
{
	int i;
	
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->Txx[i] = 0.0;
		var->Tyy[i] = 0.0;
		var->Tzz[i] = 0.0;
		
		var->Txy[i] = 0.0;
		var->Tyz[i] = 0.0;
		var->Tzx[i] = 0.0;
	}
	
	
	for (i=0; i<=local_size_ct-1; ++i){
		inter->Txx_ct[i] = 0.0;
		inter->Tyy_ct[i] = 0.0;
		inter->Tzz_ct[i] = 0.0;
	}
	
	
	for (i=0; i<=local_size_cn-1; ++i){
		inter->Txy_cn[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ued-1; ++i){
		inter->Tzx_ued[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ved-1; ++i){
		inter->Tyz_ved[i] = 0.0;
	}
	
	
}

/* ******************************************************************* */

/* ******************************************************************* */
void vel_solvers_transpose(
                           double *u,
                           double *v,
                           double *w,
                           double *rhs_u,
                           double *rhs_v,
                           double *rhs_w,
                           int substep,
                           parameters *par)
{
	
	int i, j;
	const int my_rank = par->my_rank;
	const int ny = par->ny;
	const int local_size_wp = par->local_size_wp;
	const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	
	
    
	if (par->fd_order == 2){
		vel_trisolver_v(v, rhs_v, substep, par);
		vel_trisolver_u(u, rhs_u, substep, par);
	}
	
	if (par->fd_order == 4){
		vel_pentasolver_v(v, rhs_v, substep, par);
		vel_pentasolver_u(u, rhs_u, substep, par);
	}
    
	
    
    
	if (par->const_mass == 0){
		// const pressure gradient case
		if (par->fd_order == 2){
			vel_trisolver_w(w, rhs_w, substep, par);
		}
		
		if (par->fd_order == 4){
			vel_pentasolver_w(w, rhs_w, substep, par);
		}
        
        
	}else{
		// constant mass flux case
		double *w1, *w2;
		double *rhs;
		double alpha;
		double local_bulk_w1, local_bulk_w2;
		double bulk_w1, bulk_w2;
		
		rhs = dvector(0, local_size_wp-1); // freed
		w1 = dvector(0, local_size_wp-1); // freed
		w2 = dvector(0, local_size_wp-1); // freed
		
		for (i=0; i<=local_size_wp-1; ++i){
			rhs[i] = 0.0;
			w1[i] = 0.0;
			w2[i] = 0.0;
		}
		
		
		if (par->fd_order == 2){
			vel_trisolver_w(w1, rhs_w, substep, par);
		}
		
		if (par->fd_order == 4){
			vel_pentasolver_w(w1, rhs_w, substep, par);
		}
		
		
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				rhs[par->index_wp[i][j][0]] = 1.0;
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		if (par->fd_order == 2){
			vel_trisolver_w(w2, rhs, substep, par);
		}
		if (par->fd_order == 4){
			vel_pentasolver_w(w2, rhs, substep, par);
		}
		MPI_Bcast(w2, local_size_wp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		// compute mass flux (bulk_u) with u = u1 + u2
		bulk_w1 = 0.0;
		bulk_w2 = 0.0;
		local_bulk_w1 = 0.0;
		local_bulk_w2 = 0.0;
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				local_bulk_w1 += w1[par->index_wp[i][j][0]]/(double)(local_nx*ny);
				local_bulk_w2 += w2[par->index_wp[i][j][0]]/(double)(local_nx*ny);
			}
		}
		
		MPI_Reduce(&local_bulk_w1, &bulk_w1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&local_bulk_w2, &bulk_w2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (my_rank == 0){
			bulk_w1 = bulk_w1/(double)par->nproc;
			bulk_w2 = bulk_w2/(double)par->nproc;
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Bcast(&bulk_w1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&bulk_w2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		// compute alpha for w = w1 + alpha*w2 to satisfy const par->mass_flux
		alpha = (par->mass_flux - bulk_w1)/bulk_w2;
		
		// compute w = w1 + alpha*w2
		for (i=0; i<=local_size_wp-1; ++i){
			w[i] = w1[i] + alpha*w2[i];
			rhs_w[i] += 1.0*alpha;
		}
		
		free_dvector(rhs, 0, local_size_wp-1);
		free_dvector(w1, 0, local_size_wp-1);
		free_dvector(w2, 0, local_size_wp-1);
		
    }
    
	return;
}


// getting i=end value from previous process
extern void fill_tb(double *bc, parameters *par, int tag)
{
	
	int k;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int my_rank = par->my_rank;
	
	
	MPI_Status stat1, stat2;
	double *out_r, *out_l, *in_l, *in_r;
	out_r = dvector(0, nz-1);
	out_l = dvector(0, nz-1);
	in_l  = dvector(0, nz-1);
	in_r  = dvector(0, nz-1);
	
	for (k=0; k<=nz-1; ++k){
		in_l[k]  = 0.0;
		in_r[k]  = 0.0;
		out_l[k] = 0.0;
		out_r[k] = 0.0;
	}
	
	if (my_rank == 0){
		
        for (k=0; k<=nz-1; ++k){
			out_l[k] = bc[par->index_tb[start][k]];
		}
		
		MPI_Sendrecv(out_r, nz, MPI_DOUBLE, my_rank+1, tag, \
					 in_r,  nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
        
		MPI_Sendrecv(out_l, nz, MPI_DOUBLE, par->nproc-1, tag, \
					 in_l,  nz, MPI_DOUBLE, par->nproc-1, tag, MPI_COMM_WORLD, &stat2);
        
		
	}else if (my_rank == par->nproc-1){
		
		for (k=0; k<=nz-1; ++k){
			out_l[k] = bc[par->index_tb[start][k]];
		}
		
		MPI_Sendrecv(out_l, nz, MPI_DOUBLE, my_rank-1, tag, \
					 in_l,  nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat2);
        
        MPI_Sendrecv(out_r, nz, MPI_DOUBLE, 0, tag, \
					 in_r,  nz, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat2);
		
		
	}else{
		for (k=0; k<=nz-1; ++k){
			out_l[k] = bc[par->index_tb[start][k]];
		}
		
		
		MPI_Sendrecv(out_r, nz, MPI_DOUBLE, my_rank+1, tag, \
					 in_r,  nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat1);
		
		MPI_Sendrecv(out_l, nz, MPI_DOUBLE, my_rank-1, tag, \
					 in_l,  nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat2);
		
	}
	
	
	for (k=0; k<=nz-1; ++k){
		bc[par->index_tb[end][k]] = in_r[k];
	}
	
	
	free_dvector(out_r, 0, nz-1);
	free_dvector(out_l, 0, nz-1);
	free_dvector(in_l,  0, nz-1);
	free_dvector(in_r,  0, nz-1);
	
	return;
}



