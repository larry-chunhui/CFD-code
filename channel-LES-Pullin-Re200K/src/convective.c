/*
 *  convective.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * The skew-symmetric form of the convective terms is employed.
 * Direct the simulation to appropriate subroutines including
 * the functions defined in convective_adv/div__2/4.c files.
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "share.h"
#include "interpolation.h"
#include "interpolation__2.h"
#include "interpolation__4.h"
#include "convective_div__2.h"
#include "convective_adv__2.h"
#include "convective_div__4.h"
#include "convective_adv__4.h"
#include "convective.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


void trancate_exp(double *product, char *kinds, parameters *par)
{
	
	int i, j, k;
	int index = 0;
	int x_i, x_f;
	int y_i, y_f;
	
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double m = 36.0;
	const double alpha = 0.6931*pow(1.5, m);
    double ExpDealiase[par->nz];
    
    
    
	x_i = start;
	x_f = end-1;
	y_i = 0;
	y_f = ny-1;
	
	if (strcmp(kinds, "v")==0){
		y_i = 1;
	}
	
	if (strcmp(kinds, "cn4")==0){
		x_i = start-1;
		x_f = end+1;
		y_i = -1;
		y_f = ny+1;
	}
	if (strcmp(kinds, "cn2")==0){
		x_i = start;
		x_f = end;
		y_i = 0;
		y_f = ny;
	}
	
	const double kzmax = 2*PI*(double)(nz/2)/par->Lz;
    
    
    for (k=0; k<=nz-1; ++k){
        ExpDealiase[k] = exp(-alpha*pow( par->kz[k]/kzmax, m));
    }
    
	for (i=x_i; i<=x_f; ++i){
		for (j=y_i; j<=y_f; ++j){
			for (k=0; k<=nz-1; ++k){
				
				if (strcmp(kinds, "u")==0){
					index = par->index_u[i][j][k];
				}
				if (strcmp(kinds, "v")==0){
					index = par->index_v[i][j][k];
				}
				if (strcmp(kinds, "wp")==0){
					index = par->index_wp[i][j][k];
				}
				if (strcmp(kinds, "cn2")==0 || strcmp(kinds, "cn4")==0){
					index = par->index_cn[i][j][k];
				}
				
				
                product[index] *= ExpDealiase[k];
			}
			
			if (strcmp(kinds, "u")==0){
				index = par->index_u[i][j][nz/2];
			}
			if (strcmp(kinds, "v")==0){
				index = par->index_v[i][j][nz/2];
			}
			if (strcmp(kinds, "wp")==0){
				index = par->index_wp[i][j][nz/2];
			}
			
			if (strcmp(kinds, "cn2")==0 || strcmp(kinds, "cn4")==0){
				index = par->index_cn[i][j][nz/2];
			}
			
			product[index] = 0.0;
		}
	}
	
	return;
}


void bc_trancate_exp(double *product, char *kinds, parameters *par)
{
    
	int i, k;
	int index;
	
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double m = 36.0;
	const double alpha = 0.6931*pow(1.5, m);
	const double kzmax = 2*PI*(double)(nz/2)/par->Lz;
    double ExpDealiase[par->nz];
    
    for (k=0; k<=nz-1; ++k){
        ExpDealiase[k] = exp(-alpha*pow( par->kz[k]/kzmax, m));
    }
    
    
	if (strcmp(kinds, "tb")==0 && par->les == 1){
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
                
                product[index] *= ExpDealiase[k];
			}
			index = par->index_tb[i][nz/2];
			product[index] = 0.0;
		}
	}else{
		printf("bc_trancate_exp() error");
	}
	
	return;
}








/* ******************************************************************* */
/* convective (for rhs): given u get Nu(u)							   */
/* u = vector(0, nq_x-1)											   */
/* Nu = vector(0, nq_x-1)											   */
/* v = vector(0, nq_y-1)											   */
/* Nv = vector(0, nq_y-1)											   */
/* w = vector(0, nq_w-1)											   */
/* Nw = vector(0, nq_w-1)											   */
/* ******************************************************************* */
void get_convective(
                    variables *var,
                    neighbors *shared,
                    interpolated *inter,
                    parameters *par,
                    fftwplans *planptr)
{
	
	int i, k;
	int my_rank; my_rank = par->my_rank;
	int local_nx, ny, nz;
	const int local_size_u = par->local_size_u;
	const int local_size_v = par->local_size_v;
	const int local_size_wp = par->local_size_wp;
    
    
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	
	local_nx = par->local_nx;
	ny = par->ny;
	nz = par->nz;
	
	
	double *Nu_copy, *Nv_copy, *Nw_copy;
	
	Nu_copy = FFTW_MALLOC(local_size_u);
	Nv_copy = FFTW_MALLOC(local_size_v);
	Nw_copy = FFTW_MALLOC(local_size_wp);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] = 0.0;
	}
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] = 0.0;
	}
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] = 0.0;
	}
	
	
	if (par->fd_order == 2){
		get_convective_div__2(inter, Nu_copy, Nv_copy, Nw_copy, planptr, par);
	}
	
	if (par->fd_order == 4){
		get_convective_div__4(inter, Nu_copy, Nv_copy, Nw_copy, planptr, par);
	}
    
	
	for (i=0; i<=local_size_u-1; ++i){
		var->Nu[i] = 0.5*Nu_copy[i];
	}
	for (i=0; i<=local_size_v-1; ++i){
		var->Nv[i] = 0.5*Nv_copy[i];
	}
	for (i=0; i<=local_size_wp-1; ++i){
		var->Nw[i] = 0.5*Nw_copy[i];
	}
	
	
	
	if (par->fd_order == 2){
		get_convective_adv__2(inter, Nu_copy, Nv_copy, Nw_copy, planptr, par);
	}
	
	if (par->fd_order == 4){
		get_convective_adv__4(inter, Nu_copy, Nv_copy, Nw_copy, planptr, par);
	}
	
	
	for (i=0; i<=local_size_u-1; ++i){
		var->Nu[i] += 0.5*Nu_copy[i];
	}
	for (i=0; i<=local_size_v-1; ++i){
		var->Nv[i] += 0.5*Nv_copy[i];
	}
	for (i=0; i<=local_size_wp-1; ++i){
		var->Nw[i] += 0.5*Nw_copy[i];
	}
	
	
	trancate_exp(var->Nu, "u", par);
	trancate_exp(var->Nv, "v", par);
	trancate_exp(var->Nw, "wp", par);
	
    
	fftw_free(Nu_copy);
	fftw_free(Nv_copy);
	fftw_free(Nw_copy);
	
	
	// *****************************************************************
	// for LES wall ODE
	if (par->les == 1 && par->bc == 1){
		
		bc_trancate_exp(inter->bc_duudx_t, "tb", par);
		bc_trancate_exp(inter->bc_duwdz_t, "tb", par);
		bc_trancate_exp(inter->bc_duudx_b, "tb", par);
		bc_trancate_exp(inter->bc_duwdz_b, "tb", par);
		
		bc_trancate_exp(inter->bc_ududx_t, "tb", par);
		bc_trancate_exp(inter->bc_wdudz_t, "tb", par);
		bc_trancate_exp(inter->bc_ududx_b, "tb", par);
		bc_trancate_exp(inter->bc_wdudz_b, "tb", par);
		
		
		int index;
		int index1, index2, index3, index4;
		const int j_log = par->j_log;
		
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
				index1 = par->index_les[i][j_log][k];
				index2 = par->index_les[i][ny-j_log][k];
				index3 = par->index_les[i][j_log+1][k];
				index4 = par->index_les[i][ny-1-j_log][k];
				
				if (par->ode_choice == 1) {
					inter->bc_uv_b[index]
					= (inter->u_phys_les[index1]*inter->v_phys_les[index1]
					   + inter->u_phys_les[index3]*inter->v_phys_les[index3])/2.0;
                    
					inter->bc_uv_t[index]
					= (inter->u_phys_les[index2]*inter->v_phys_les[index2]
					   + inter->u_phys_les[index4]*inter->v_phys_les[index4])/2.0;
				}
				
				if (par->ode_choice == 2) {
					inter->bc_uv_b[index]
					= inter->u_phys_les[index3]*inter->v_phys_les[index3];
					
					inter->bc_uv_t[index]
					= inter->u_phys_les[index4]*inter->v_phys_les[index4];
				}
				
			}
		}
		
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_uv_t, inter->bc_uv_t);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_uv_b, inter->bc_uv_b);
		
		for (i=0; i<=(local_nx+1)*nz-1; ++i){
			inter->bc_uv_b[i] = inter->bc_uv_b[i]/(double)nz; // scale
			inter->bc_uv_t[i] = inter->bc_uv_t[i]/(double)nz; // scale
		}
		
		bc_trancate_exp(inter->bc_uv_t, "tb", par);
		bc_trancate_exp(inter->bc_uv_b, "tb", par);
		
	}
    
	
	return;
}






/* ******************************************************************* */
/*
 product uu, vv and ww defined at center of each cell
 uu[]; index_ct[][][];[start-1->end][-1->ny][0->nz-1]
 */
/* ******************************************************************* */
void calc_multiply_ct(
                      double *u_ct1,
                      double *u_ct2,
                      double *uu,
                      parameters *par)
{
	int i;
	const int local_size_ct = par->local_size_ct;
    
	for (i=0; i<=local_size_ct-1; ++i){
		uu[i] = u_ct1[i]*u_ct2[i];
	}
	
	return;
}


/* ******************************************************************* */
/*
 product uw defined at edges ued[]
 uw[]; index_ued[][][] [start-1->end][-1->ny][0->nz-1]
 */
/* ******************************************************************* */
void calc_multiply_ued(
                       double *u_ed,
                       double *w_ued,
                       double *uw,
                       parameters *par)
{
	int i;
	const int local_size_ued = par->local_size_ued;
    
	for (i=0; i<=local_size_ued-1; ++i){
		uw[i] = u_ed[i]*w_ued[i];
	}
    
	return;
}


/* ******************************************************************* */
/*
 product vw defined at edges ved[]
 vw[]; index_ved[][][] [start-1->end][0->ny][0->nz-1]
 */
/* ******************************************************************* */
void calc_multiply_ved(
                       double *v_ed,
                       double *w_ved,
                       double *vw,
                       parameters *par)
{
	int i;
	const int local_size_ved = par->local_size_ved;
    
	for (i=0; i<=local_size_ved-1; ++i){
		vw[i] = v_ed[i]*w_ved[i];
	}
	
	return;
}


/* ******************************************************************* */
/*
 product uv defined at corners
 uv[]; [start->end][0->ny][0->nz-1] index_cn[][][]
 */
/* ******************************************************************* */
void calc_multiply_cn(
                      double *u_cn,
                      double *v_cn,
                      double *uv,
                      parameters *par)
{
	int i;
	const int local_size_cn = par->local_size_cn;
    
	for (i=0; i<=local_size_cn-1; ++i){
		uv[i] = u_cn[i]*v_cn[i];
	}
    
	return;
}

