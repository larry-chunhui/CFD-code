/*
 *  near_wall_ode.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Soves the ODE for dudy -- first part of the wall model.
 * Details found is Chung & Pullin 2009, Inoue and Pullin 2012 and Saito et al. 2012.
 * Roughness correction with Colebrook's formula -- Saito et al. 2012.
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "definitions.h"
#include "nrutil.h"
#include "convective.h"
#include "differentiate.h"
#include "near_wall_ode.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


extern void bc__time_march_bc_dudy(
								   variables *var_prev,
								   variables *var_curr,
								   variables *var_next,
								   parameters *par )
{
	int i, k;
	int index;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
			
			var_next->bc_dudy_t[index]
			= var_curr->bc_rhs_dudy_t[index];
			
			var_next->bc_dudy_b[index]
			= var_curr->bc_rhs_dudy_b[index];
			
		}
	}
	
	
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
			
			var_prev->Hbc_dudy_t[index]
			= var_curr->Hbc_dudy_t[index];
			
			var_prev->Hbc_dudy_b[index]
			= var_curr->Hbc_dudy_b[index];
		}
	}
	
	return;
}



extern void bc__time_march_bc_dudy_analitic(
                                            variables *var0,
                                            variables *var1,
                                            int substep,
                                            parameters *par )
{
	
	int i;
	const int local_nx = par->local_nx;
	const int nz = par->nz;
	
	const double gamma = par->rk3_gamma[substep];
	const double zeta = par->rk3_zeta[substep];
	const double dt_substep = (gamma+zeta)*par->dt;
	
	
	double scratch;
	
	for (i=0; i<=(local_nx+1)*(nz)-1; ++i){
		scratch = 1.0/var0->eta0_b[i]*          exp( - dt_substep*var0->lambda_b[i]*var0->eta0_bar_b[i] )
		+ 1.0/var0->eta0_bar_b[i]*(1.0- exp( - dt_substep*var0->lambda_b[i]*var0->eta0_bar_b[i] ));
		var1->bc_dudy_b[i] = 1.0/scratch;
	}
	
	
	for (i=0; i<=(local_nx+1)*(nz)-1; ++i){
		scratch = 1.0/var0->eta0_t[i]*          exp( - dt_substep*var0->lambda_t[i]*var0->eta0_bar_t[i] )
		+ 1.0/var0->eta0_bar_t[i]*(1.0- exp( - dt_substep*var0->lambda_t[i]*var0->eta0_bar_t[i] ));
		var1->bc_dudy_t[i] = 1.0/scratch;
	}
    
	
	
}



extern void term_check__nonlinear_ode(
									  variables *var,
									  interpolated *inter,
									  int substep,
									  parameters *par,
									  fftwplans *planptr)
{
	int i, k;
	int index;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int j_log = par->j_log;
	
	
	double *term_A1;
	double *term_A2;
	double *term_B1;
	double *term_B2;
	double *term_C1;
	double *term_C2;
	double *term_D;
	double *term_E1;
	double *term_E2;
	
	term_A1 = FFTW_MALLOC(local_nx);
	term_A2 = FFTW_MALLOC(local_nx);
	term_B1 = FFTW_MALLOC(local_nx);
	term_B2 = FFTW_MALLOC(local_nx);
	term_C1 = FFTW_MALLOC(local_nx);
	term_C2 = FFTW_MALLOC(local_nx);
	term_D  = FFTW_MALLOC(local_nx);
	term_E1 = FFTW_MALLOC(local_nx);
	term_E2 = FFTW_MALLOC(local_nx);
	
	
    double bc_h_log[2] = {0.0, 0.0};  // Distance from log plane to physical wall.
	
	if (par->ode_choice == 1) {
		bc_h_log[0] = ( par->h0 + ((double)j_log+0.5)*par->dy );
		bc_h_log[1] = - ( par->h0 + ((double)j_log+0.5)*par->dy );
	}
	
	if (par->ode_choice == 2) {
		bc_h_log[0] = ( par->h0 + ((double)j_log+1.0)*par->dy );
		bc_h_log[1] = - ( par->h0 + ((double)j_log+1.0)*par->dy );
	}
	
	double *scratch3;
	double *scratch4;
	
	scratch3 = FFTW_MALLOC(par->local_size_u);
    
    if (par->fd_order == 2){
		derivative_x1_ct2u(inter->Txx_ct, scratch3, par);
	}
    
    if (par->fd_order == 4){
        scratch4 = FFTW_MALLOC(par->local_size_u);
        
        derivative_x1_ct2u(inter->Txx_ct, scratch3, par);
        derivative_x3_ct2u(inter->Txx_ct, scratch4, par);
        
        // data in Fourier space
        for (i=start; i<=end-1; ++i) {
            for (k=0; k<=nz-1; ++k) {
                index = par->index_tb[i][k];
                if (par->ode_choice == 1) {
                    
                    // Skew-symmetric form
                    scratch3[par->index_u[i][j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][j_log][k]];
                }
                
                if (par->ode_choice == 2) {
                    
                    // Skew-symmetric form
                    scratch3[par->index_u[i][j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][j_log][k]];
                    
                    scratch3[par->index_u[i][j_log-1][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][j_log+1][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][j_log+1][k]];
                    
                }
            }
        }
        fftw_free(scratch4);
    }
    
    
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		index = par->index_tb[i][0];
		
		if (par->ode_choice == 1) {
			// Skew-symmetric form
			term_A1[i-start] = (inter->bc_duudx_b[index] + inter->bc_ududx_b[index])/2.0;
			term_A2[i-start] = scratch3[par->index_u[i][j_log][0]];
		}
		if (par->ode_choice == 2) {
			
			// Skew-symmetric form
			term_A1[i-start] = (inter->bc_duudx_b[index] + inter->bc_ududx_b[index])/2.0;
			term_A2[i-start] = (scratch3[par->index_u[i][j_log][0]] + scratch3[par->index_u[i][j_log+1][0]])/2.0;
		}
	}
	
	
	
	
	
	derivative_z_ued2u(inter->Tzx_ued, scratch3, par);
	
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		index = par->index_tb[i][0];
		
		if (par->ode_choice == 1) {
			
			// Skew-symmetric form
			term_B1[i-start] = (inter->bc_duwdz_b[index] + inter->bc_wdudz_b[index])/2.0;
			term_B2[i-start] = scratch3[par->index_u[i][j_log][0]];
		}
		
		
		if (par->ode_choice == 2) {
			
			// Skew-symmetric form
			term_B1[i-start] = (inter->bc_duwdz_b[index] + inter->bc_wdudz_b[index])/2.0;
			term_B2[i-start] = (scratch3[par->index_u[i][j_log][0]] + scratch3[par->index_u[i][j_log+1][0]])/2.0;
		}
	}
	
	
	fftw_free(scratch3);
	
	
	double h_log = bc_h_log[0];
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		
		index = par->index_tb[i][0];
		
		if (par->ode_choice == 1) {
			
			term_C1[i-start] = - inter->bc_uv_b[index]/h_log;
			term_C2[i-start] = - (inter->Txy_cn[par->index_cn[i][j_log][0]]
								  + inter->Txy_cn[par->index_cn[i][j_log+1][0]])/2.0/h_log;
		}
		
		if (par->ode_choice == 2) {
			
			term_C1[i-start] = - inter->bc_uv_b[index]/h_log;
			term_C2[i-start] = - inter->Txy_cn[par->index_cn[i][j_log+1][0]]/h_log;
			
		}
		
	}
	
	
	
	for (i=start; i<=end-1; ++i) {
		// bottom wall
		if (par->ode_choice == 1) {
			term_D[i-start] = - var->Gpx[par->index_u[i][j_log][0]] - par->dPdx;
		}
		if (par->ode_choice == 2) {
			term_D[i-start] = - (var->Gpx[par->index_u[i][j_log][0]]+var->Gpx[par->index_u[i][j_log+1][0]])/2.0 - par->dPdx;
		}
	}
	
	
	
	
	
	double *scratch_b;
	scratch_b = FFTW_MALLOC((local_nx+1)*nz);
	
	int ijk_u, ijk_cn;
	const double nu = 1.0/par->re;
	
	for (i=0; i<=(local_nx+1)*nz-1; ++i){
		scratch_b[i] = var->bc_dudy_b[i]; // physical space
	}
	
	fftw_execute_r2r(planptr->p1d_z_bc_tb, scratch_b, scratch_b);
	
	for (i=0; i<=(local_nx+1)*nz-1; ++i){
		scratch_b[i] = scratch_b[i]/(double)nz; //rescale
	}
	
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		// bottom boundary
		index = par->index_tb[i][0];
		
		term_E2[i-start] = - nu/h_log*scratch_b[index];
	}
	
	
	
	// physical space
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			// bottom wall
			ijk_u = par->index_u[i][j_log][k];
			ijk_cn = par->index_cn[i][j_log+1][k];
			index = par->index_tb[i][k];
			
			if (par->ode_choice == 1) {
				scratch_b[index] = (var->dudy[par->index_cn[i][j_log+1][k]] + var->dudy[par->index_cn[i][j_log][k]])/2.0;
			}
			if (par->ode_choice == 2) {
				scratch_b[index] = var->dudy[ijk_cn];
			}
		}
	}
	
	fftw_execute_r2r(planptr->p1d_z_bc_tb, scratch_b, scratch_b);
	
	for (i=0; i<=(local_nx+1)*nz-1; ++i){
		scratch_b[i] = scratch_b[i]/(double)nz; //rescale
	}
	
	
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		// bottom boundary
		index = par->index_tb[i][0];
		
		term_E1[i-start] = nu/h_log*scratch_b[index];
	}
	
	
	fftw_free(scratch_b);
	
	par->A1 = average__nonlinear_ode(term_A1, par);
	par->A2 = average__nonlinear_ode(term_A2, par);
	par->B1 = average__nonlinear_ode(term_B1, par);
	par->B2 = average__nonlinear_ode(term_B2, par);
	par->C1 = average__nonlinear_ode(term_C1, par);
	par->C2 = average__nonlinear_ode(term_C2, par);
	par->D  = average__nonlinear_ode(term_D, par);
	par->E1 = average__nonlinear_ode(term_E1, par);
	par->E2 = average__nonlinear_ode(term_E2, par);
	
	
	mean_values_1d__term_check_nonlinear(term_A1, var->A1_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_A2, var->A2_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_B1, var->B1_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_B2, var->B2_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_C1, var->C1_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_C2, var->C2_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_D,  var->D_mean,  substep, par);
	mean_values_1d__term_check_nonlinear(term_E1, var->E1_mean, substep, par);
	mean_values_1d__term_check_nonlinear(term_E2, var->E2_mean, substep, par);
	
	{
		
		if (par->it%par->stat == 0 && substep == 0){
			wall_terms_check_dump(term_A1, term_A2,
								  term_B1, term_B2,
								  term_C1, term_C2,
								  term_D,
								  term_E1, term_E2,
								  "snap", par);
			wall_terms_check_dump(var->A1_mean, var->A2_mean,
								  var->B1_mean, var->B2_mean,
								  var->C1_mean, var->C2_mean,
								  var->D_mean,
								  var->E1_mean, var->E2_mean,
								  "mean", par);
		}
		
		
		
	}
	
	
	fftw_free(term_A1);
	fftw_free(term_A2);
	fftw_free(term_B1);
	fftw_free(term_B2);
	fftw_free(term_C1);
	fftw_free(term_C2);
	fftw_free(term_D);
	fftw_free(term_E1);
	fftw_free(term_E2);
	
	
	return;
}


extern void mean_values_1d__term_check_nonlinear(
												 double *snap,
												 double *mean,
												 int substep,
												 parameters *par)
{
	
	int i;
	
	if (par->it == 0 && substep == 0){
		for (i=0; i<=par->local_nx-1; ++i){
			mean[i] = snap[i];
		}
	}
	
	
	double copy;
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
	
	const double dt = (alpha+beta)*par->dt;
	const double T_ave = 2.0;
	
	for (i=0; i<=par->local_nx-1; ++i){
		copy = mean[i];
		mean[i] = dt/T_ave*snap[i] + (1.0 - dt/T_ave)*copy;
	}
	
	return;
}


extern void wall_terms_check_dump(
								  double *A1, double *A2,
								  double *B1, double *B2,
								  double *C1, double *C2,
								  double *D,
								  double *E1, double *E2,
								  const char *basename,
								  parameters *par)
{
	
	
	int nx, ny, nz, total_local_size, total_size;
    double *a1, *a2, *b1, *b2, *c1, *c2, *d, *e1, *e2;
	
	nx = par->nx;
    ny = par->ny;
    nz = par->nz;
	
    total_local_size = par->local_nx;
    total_size       = par->local_nx*par->nproc;
	
    if ( par->my_rank == 0 ) {
        a1 = ( double * )fftw_malloc( sizeof( double )*total_size );
		a2 = ( double * )fftw_malloc( sizeof( double )*total_size );
		b1 = ( double * )fftw_malloc( sizeof( double )*total_size );
		b2 = ( double * )fftw_malloc( sizeof( double )*total_size );
		c1 = ( double * )fftw_malloc( sizeof( double )*total_size );
		c2 = ( double * )fftw_malloc( sizeof( double )*total_size );
		d = ( double * )fftw_malloc( sizeof( double )*total_size );
		e1 = ( double * )fftw_malloc( sizeof( double )*total_size );
		e2 = ( double * )fftw_malloc( sizeof( double )*total_size );
    }
    else {
        a1 = NULL;
		a2 = NULL;
		b1 = NULL;
		b2 = NULL;
		c1 = NULL;
		c2 = NULL;
		d = NULL;
		e1 = NULL;
		e2 = NULL;
    }
	
    MPI_Gatherv(A1,     total_local_size,         MPI_DOUBLE,
				a1,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(A2,     total_local_size,         MPI_DOUBLE,
				a2,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(B1,     total_local_size,         MPI_DOUBLE,
				b1,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(B2,     total_local_size,         MPI_DOUBLE,
				b2,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(C1,     total_local_size,         MPI_DOUBLE,
				c1,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(C2,     total_local_size,         MPI_DOUBLE,
				c2,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(D,     total_local_size,         MPI_DOUBLE,
				d,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(E1,     total_local_size,         MPI_DOUBLE,
				e1,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	MPI_Gatherv(E2,     total_local_size,         MPI_DOUBLE,
				e2,  par->recvcounts_mean_tb, par->displs_mean_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	
    if ( par->my_rank == 0 ) {
		
        FILE *fp;
        char filename[128];
        int i;
		double x;
		
        sprintf( filename, "./outputdir/wall_terms_check_%s_it%d.txt",
				basename,
				par->it );
        fp = fopen( filename, "w" );
        if ( fp == NULL ) {
			printf( "wall_terms_check_dump(): cannot write to file" );
		}
		
		for (i=0; i<=(par->local_nx)*par->nproc-1; ++i){
			x = i*par->dx;
			fprintf( fp, "%+e	%+e	%+e	%+e	%+e	%+e	%+e	%+e	%+e	%+e\n",
					x, a1[i], a2[i], b1[i], b2[i], c1[i], c2[i], d[i], e1[i], e2[i]);
		}
		
        fclose( fp );
    }
	
    if ( par->my_rank == 0 ) {
        fftw_free( a1 );
		fftw_free( a2 );
		fftw_free( b1 );
		fftw_free( b2 );
		fftw_free( c1 );
		fftw_free( c2 );
		fftw_free( d );
		fftw_free( e1 );
		fftw_free( e2 );
    }
	
	return;
}



double average__nonlinear_ode(double *term, parameters *par)
{
	int i;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	// check the average of each terms
	
	double sum;
	double local_sum;
	double average;
	local_sum = 0.0;
	
	for (i=start; i<=end-1; ++i) {
		local_sum += term[i-start];
	}
	
	MPI_Reduce(&local_sum, &sum,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	if (par->my_rank == 0) {
		average = sum/par->nx;
	}
	
	MPI_Bcast( &average,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	
	
	return average;
}



extern void bc__nonlinear(
						  variables *var,
						  interpolated *inter,
						  int substep,
						  parameters *par,
						  fftwplans *planptr)
{
    
	int i, k;
	int index;
	const int nz = par->nz;
	const int ny = par->ny;
	const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int j_log = par->j_log;
	
	double *scratch_t, *scratch_b;
	scratch_t = FFTW_MALLOC((local_nx+1)*nz);
	scratch_b = FFTW_MALLOC((local_nx+1)*nz);
	
    double bc_h_log[2] = {0.0, 0.0};  // Distance from log plane to physical wall.
    
	
	if (par->ode_choice == 1) {
		bc_h_log[0] = ( par->h0 + ((double)j_log+0.5)*par->dy );
		bc_h_log[1] = - ( par->h0 + ((double)j_log+0.5)*par->dy );
	}
	
	if (par->ode_choice == 2) {
		bc_h_log[0] = ( par->h0 + ((double)j_log+1.0)*par->dy );
		bc_h_log[1] = - ( par->h0 + ((double)j_log+1.0)*par->dy );
	}
	
	double *scratch3;
	double *scratch4;
	
	scratch3 = FFTW_MALLOC(par->local_size_u);
    
    if (par->fd_order == 2){
		derivative_x1_ct2u(inter->Txx_ct, scratch3, par);
	}
    
    if (par->fd_order == 4){
        scratch4 = FFTW_MALLOC(par->local_size_u);
        
        derivative_x1_ct2u(inter->Txx_ct, scratch3, par);
        derivative_x3_ct2u(inter->Txx_ct, scratch4, par);
        
        // data in Fourier space
        for (i=start; i<=end-1; ++i) {
            for (k=0; k<=nz-1; ++k) {
                index = par->index_tb[i][k];
                if (par->ode_choice == 1) {
                    
                    // Skew-symmetric form
                    scratch3[par->index_u[i][ny-1-j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][ny-1-j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][ny-1-j_log][k]];
                    
                    scratch3[par->index_u[i][j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][j_log][k]];
                }
                
                if (par->ode_choice == 2) {
                    
                    // Skew-symmetric form
                    scratch3[par->index_u[i][ny-1-j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][ny-1-j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][ny-1-j_log][k]];
                    
                    scratch3[par->index_u[i][j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][j_log][k]];
                    
                    scratch3[par->index_u[i][ny-2-j_log][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][ny-2-j_log][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][ny-2-j_log][k]];
                    
                    scratch3[par->index_u[i][j_log-1][k]]
                    = 9.0/8.0*scratch3[par->index_u[i][j_log+1][k]]
                    - 1.0/8.0*scratch4[par->index_u[i][j_log+1][k]];
                    
                }
            }
        }
        fftw_free(scratch4);
    }
	
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
			if (par->ode_choice == 1) {
				
				// Skew-symmetric form
				var->bc_convective_t[index]
				= (inter->bc_duudx_t[index] + inter->bc_ududx_t[index])/2.0
				+ scratch3[par->index_u[i][ny-1-j_log][k]];
				
				var->bc_convective_b[index]
				= (inter->bc_duudx_b[index] + inter->bc_ududx_b[index])/2.0
				+ scratch3[par->index_u[i][j_log][k]];
				
				/*
                 // divergence form
                 var->bc_convective_t[index]
                 = inter->bc_duudx_t[index]
                 + scratch3[par->index_u[i][ny-1-j_log][k]];
                 
                 var->bc_convective_b[index]
                 = inter->bc_duudx_b[index]
                 + scratch3[par->index_u[i][j_log][k]];
                 */
			}
			
			if (par->ode_choice == 2) {
				
				// Skew-symmetric form
				var->bc_convective_t[index]
				= (inter->bc_duudx_t[index] + inter->bc_ududx_t[index])/2.0
				+ (scratch3[par->index_u[i][ny-1-j_log][k]] + scratch3[par->index_u[i][ny-2-j_log][k]])/2.0;
				
				var->bc_convective_b[index]
				= (inter->bc_duudx_b[index] + inter->bc_ududx_b[index])/2.0
				+ (scratch3[par->index_u[i][j_log][k]] + scratch3[par->index_u[i][j_log+1][k]])/2.0;
				
				/*
                 // divergence form
                 var->bc_convective_t[index]
                 = inter->bc_duudx_t[index]
                 + (scratch3[par->index_u[i][ny-1-j_log][k]] + scratch3[par->index_u[i][ny-2-j_log][k]])/2.0;
                 
                 var->bc_convective_b[index]
                 = inter->bc_duudx_b[index]
                 + (scratch3[par->index_u[i][j_log][k]] + scratch3[par->index_u[i][j_log+1][k]])/2.0;
                 */
			}
		}
	}
	
	
	derivative_z_ued2u(inter->Tzx_ued, scratch3, par);
	
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
			if (par->ode_choice == 1) {
				
				// Skew-symmetric form
				var->bc_convective_t[index]
				+= (inter->bc_duwdz_t[index] + inter->bc_wdudz_t[index])/2.0
				+ scratch3[par->index_u[i][ny-1-j_log][k]];
				
				
				var->bc_convective_b[index]
				+= (inter->bc_duwdz_b[index] + inter->bc_wdudz_b[index])/2.0
				+ scratch3[par->index_u[i][j_log][k]];
				
				/*
                 // divergence form
                 var->bc_convective_t[index]
                 += inter->bc_duwdz_t[index]
                 + scratch3[par->index_u[i][ny-1-j_log][k]];
                 
                 var->bc_convective_b[index]
                 += inter->bc_duwdz_b[index]
                 + scratch3[par->index_u[i][j_log][k]];
                 */
			}
			
			
			if (par->ode_choice == 2) {
				
				// Skew-symmetric form
				var->bc_convective_t[index]
				+= (inter->bc_duwdz_t[index] + inter->bc_wdudz_t[index])/2.0
				+ (scratch3[par->index_u[i][ny-1-j_log][k]] + scratch3[par->index_u[i][ny-2-j_log][k]])/2.0;
				
				var->bc_convective_b[index]
				+= (inter->bc_duwdz_b[index] + inter->bc_wdudz_b[index])/2.0
				+ (scratch3[par->index_u[i][j_log][k]] + scratch3[par->index_u[i][j_log+1][k]])/2.0;
				
				/*
                 // divergence form
                 var->bc_convective_t[index]
                 += inter->bc_duwdz_t[index]
                 + (scratch3[par->index_u[i][ny-1-j_log][k]] + scratch3[par->index_u[i][ny-2-j_log][k]])/2.0;
                 
                 
                 var->bc_convective_b[index]
                 += inter->bc_duwdz_b[index]
                 + (scratch3[par->index_u[i][j_log][k]] + scratch3[par->index_u[i][j_log+1][k]])/2.0;
                 */
			}
            
		}
	}
	
	
	fftw_free(scratch3);
	
	// check a sum of convective terms
	{
		double sum_t, sum_b;
		double local_sum_t, local_sum_b;
		local_sum_t = 0.0;
		local_sum_b = 0.0;
		
		for (i=start; i<=end-1; ++i) {
			index = par->index_tb[i][0];
			local_sum_t += var->bc_convective_t[index];
			local_sum_b += var->bc_convective_b[index];
		}
        
		MPI_Reduce(&local_sum_t, &sum_t,
				   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce(&local_sum_b, &sum_b,
				   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		
		if (par->my_rank == 0) {
			sum_t = sum_t/par->nx*(par->h0+0.5*par->dz);
			sum_b = sum_b/par->nx*(par->h0+0.5*par->dz);
			printf("at substep %d sum of comvective terms top(bot) %+e	(%+e)\n", substep, sum_t, sum_b);
		}
		
	}
	
	fftw_execute_r2r(planptr->p1d_invz_bc_tb, var->bc_convective_t, var->bc_convective_t);
	fftw_execute_r2r(planptr->p1d_invz_bc_tb, var->bc_convective_b, var->bc_convective_b);
	
	double *Gpx_copy;
	Gpx_copy = FFTW_MALLOC(par->local_size_u);
	
	for (i=0; i<=par->local_size_u-1; ++i){
		Gpx_copy[i] = var->Gpx[i];
	}
	
	fftw_execute_r2r(planptr->p1d_invz_u, Gpx_copy, Gpx_copy);
	
	
	// data in Fourier space
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			
			index = par->index_tb[i][k];
			
			if (par->ode_choice == 1) {
				
				scratch_t[index] = inter->bc_uv_t[index]
				+ (inter->Txy_cn[par->index_cn[i][ny-j_log][k]]
				   + inter->Txy_cn[par->index_cn[i][ny-1-j_log][k]])/2.0;
                
				scratch_b[index] = inter->bc_uv_b[index]
				+ (inter->Txy_cn[par->index_cn[i][j_log][k]]
				   + inter->Txy_cn[par->index_cn[i][j_log+1][k]])/2.0;
			}
			
			if (par->ode_choice == 2) {
				
				scratch_t[index] = inter->bc_uv_t[index]
				+ inter->Txy_cn[par->index_cn[i][ny-1-j_log][k]];
				
				scratch_b[index] = inter->bc_uv_b[index]
				+ inter->Txy_cn[par->index_cn[i][j_log+1][k]];
			}
            
		}
	}
	
	fftw_execute_r2r(planptr->p1d_invz_bc_tb, scratch_t, scratch_t);
	fftw_execute_r2r(planptr->p1d_invz_bc_tb, scratch_b, scratch_b);
	
	
	{
		int ijk_u, ijk_cn;
		const double nu = 1.0/par->re;
		double dpdx_log = 0.0, h_log, dudy_log = 0.0, uv_log;
		double convective;
        
		for (i=start; i<=end-1; ++i) {
			for (k=0; k<=nz-1; ++k) {
				// bottom wall
				ijk_u = par->index_u[i][j_log][k];
				ijk_cn = par->index_cn[i][j_log+1][k];
				index = par->index_tb[i][k];
                
				h_log = bc_h_log[0];
				if (par->ode_choice == 1) {
					dudy_log = (var->dudy[par->index_cn[i][j_log+1][k]] + var->dudy[par->index_cn[i][j_log][k]])/2.0;
					dpdx_log = Gpx_copy[ijk_u];
				}
				if (par->ode_choice == 2) {
					dudy_log = var->dudy[ijk_cn];
					dpdx_log = (Gpx_copy[par->index_u[i][j_log][k]]+Gpx_copy[par->index_u[i][j_log+1][k]])/2.0;
				}
				dpdx_log += par->dPdx; // sym checked
				uv_log = scratch_b[index]; // sym checked
				convective = var->bc_convective_b[index]; // sym checked
                
				var->Hbc_dudy_b[index] = -uv_log + nu*dudy_log
				- h_log*( dpdx_log + convective);
				
				var->eta0_bar_b[index] = var->Hbc_dudy_b[index]/nu;
                
				// top wall
				ijk_u = par->index_u[i][ny-1-j_log][k];
				ijk_cn = par->index_cn[i][ny-1-j_log][k];
				index = par->index_tb[i][k];
                
				h_log = bc_h_log[1];
				if (par->ode_choice == 1) {
					dudy_log = (var->dudy[par->index_cn[i][ny-j_log][k]] + var->dudy[par->index_cn[i][ny-1-j_log][k]])/2.0;
					dpdx_log = Gpx_copy[ijk_u];
				}
				if (par->ode_choice == 2) {
					dudy_log = var->dudy[ijk_cn]; // sym checked
					dpdx_log = (Gpx_copy[par->index_u[i][ny-1-j_log][k]]+Gpx_copy[par->index_u[i][ny-2-j_log][k]])/2.0;
				}
				dpdx_log += par->dPdx; // sym checked
				uv_log = scratch_t[index]; // sym checked
				convective = var->bc_convective_t[index]; // sym checked
                
				var->Hbc_dudy_t[index] = - uv_log + nu*dudy_log
				- h_log*( dpdx_log + convective);
                
				var->eta0_bar_t[index] = var->Hbc_dudy_t[index]/nu;
			}
		}
	}
    
	fftw_free(Gpx_copy);
	
	double *ode_roughness_b, *ode_roughness_t;
	ode_roughness_b = FFTW_MALLOC((par->local_nx+1)*(par->nz));
	ode_roughness_t = FFTW_MALLOC((par->local_nx+1)*(par->nz));
	bc_rhs_roughness(var, ode_roughness_b, ode_roughness_t, par);
    
	{
		int ijk_cn, ijk_les, ijk_ued;
		const double nu = 1.0/par->re;
		double bc_dudy, u_log = 0.0, h_log;
		
		for (i=0; i<=(local_nx+1)*nz-1; ++i){
			scratch_t[i] = var->bc_dudy_t[i]; // physical space
			scratch_b[i] = var->bc_dudy_b[i]; // physical space
		}
		
		
		// data in physical space
		for (i=start; i<=end-1; ++i) {
			for (k=0; k<=nz-1; ++k) {
				// bottom boundary
				ijk_les = par->index_les[i][j_log+1][k];
				ijk_ued = par->index_ued[i][j_log][k];
				ijk_cn = par->index_cn[i][j_log+1][k];
				index = par->index_tb[i][k];
				
				h_log = bc_h_log[0];
				bc_dudy = scratch_b[index];
				
				
				if (par->ode_choice == 1) {
					u_log = inter->u_phys_ued[ijk_ued];
				}
				if (par->ode_choice == 2) {
					u_log = inter->u_phys_les[ijk_les];
				}
				
				var->eta0_b[index] = scratch_b[index];
				var->lambda_b[index] = 2.0*nu/(u_log*h_log);
				
                
				var->Hbc_dudy_b[index]
                = 2.0*(var->Hbc_dudy_b[index]-nu*bc_dudy)/(u_log/bc_dudy-sqrt(nu/fabs(bc_dudy))*ode_roughness_b[index])/h_log;
				
				
				// top boundary
				ijk_les = par->index_les[i][ny-1-j_log][k];
				ijk_ued = par->index_ued[i][ny-1-j_log][k];
				ijk_cn = par->index_cn[i][ny-1-j_log][k];
				index = par->index_tb[i][k];
				
				h_log = bc_h_log[1];
				bc_dudy = scratch_t[index];
				
				
				if (par->ode_choice == 1) {
					u_log = inter->u_phys_ued[ijk_ued];
				}
				if (par->ode_choice == 2) {
					u_log = inter->u_phys_les[ijk_les];
				}
				
				var->eta0_t[index] = scratch_t[index];
				var->lambda_t[index] = 2.0*nu/(u_log*h_log);
				
                var->Hbc_dudy_t[index]
                = 2.0*(var->Hbc_dudy_t[index]-nu*bc_dudy)/(u_log/bc_dudy+sqrt(nu/fabs(bc_dudy))*ode_roughness_t[index])/h_log;
				
			}
		}
		
	}
	
	
	fftw_free(scratch_b);
	fftw_free(scratch_t);
	fftw_free(ode_roughness_t);
	fftw_free(ode_roughness_b);
	return;
}



extern void bc__rhs(
					variables *var_prev,
					variables *var_curr,
					int substep,
					parameters *par )
{
	int i, k;
	int index;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double dt = par->dt;
	const double gamma = par->rk3_gamma[substep];
	const double zeta = par->rk3_zeta[substep];
	
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
			
			var_curr->bc_rhs_dudy_t[index]
			= var_curr->bc_dudy_t[index]
			+ dt*( gamma*var_curr->Hbc_dudy_t[index]
                  + zeta*var_prev->Hbc_dudy_t[index] );
			
			var_curr->bc_rhs_dudy_b[index]
			= var_curr->bc_dudy_b[index]
			+ dt*( gamma*var_curr->Hbc_dudy_b[index]
                  + zeta*var_prev->Hbc_dudy_b[index] );
		}
    }
	
	return;
}


extern void bc_rhs_roughness(variables *var, double *ode_roughness_b, double *ode_roughness_t, parameters *par){
    int bc_size = (par->local_nx+1)*(par->nz);
    int i;
    double retau_b, retau_t;   // Dynamically calculated Re_tua
    double nu = 1.0/par->re;
    for (i=0; i<bc_size; i++){
        retau_b=sqrt(fabs(var->bc_dudy_b[i]/nu));   // = delta*u_tau_dyn/nu but u_tau = sqrt(eta*nu)
        retau_t=sqrt(fabs(var->bc_dudy_t[i]/nu));   // so = 1*sqrt(eta/nu)
        ode_roughness_b[i]=0.26/var->bc_kappa_b[i]*par->eps*retau_b/(1.0+0.26*par->eps*retau_b);
        ode_roughness_t[i]=0.26/var->bc_kappa_t[i]*par->eps*retau_t/(1.0+0.26*par->eps*retau_t);
    }
    return;
}
