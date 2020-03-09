/*
 *  convective_div__4.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Calculates the divergence form of the convective terms with
 * fourth order schemes.
 */

#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "share.h"
#include "convective.h"
#include "convective_div__4.h"
#include "differentiate.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )



/* ******************************************************************* */
/* convective (for rhs): given u get Nu(u)							   */
/* u = vector(0, nq_x-1)											   */
/* Nu = vector(0, nq_x-1)											   */
/* v = vector(0, nq_y-1)											   */
/* Nv = vector(0, nq_y-1)											   */
/* w = vector(0, nq_w-1)											   */
/* Nw = vector(0, nq_w-1)											   */
/* ******************************************************************* */
void get_convective_div__4(
                           interpolated *inter,
                           double *Nu_copy,
                           double *Nv_copy,
                           double *Nw_copy,
                           fftwplans *planptr,
                           parameters *par)
{
	
	get_convective_div__4_Nu(inter, Nu_copy, planptr, par);
	get_convective_div__4_Nv(inter, Nv_copy, planptr, par);
	get_convective_div__4_Nw(inter, Nw_copy, planptr, par);
    
	return;
}




/* ******************************************************************* */
/* convective (for rhs): given u get Nu(u)							   */
/* u = vector(0, nq_x-1)											   */
/* Nu = vector(0, nq_x-1)											   */
/* ******************************************************************* */
void get_convective_div__4_Nu(
                              interpolated *inter,
                              double *Nu_copy,
                              fftwplans *planptr,
                              parameters *par)
{
	int i, k, index;
	int index1, index2;
	int index3, index4;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int local_size_u = par->local_size_u;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	double *uu1, *uu3;
	double *uv1, *uv3;
	double *uw;
	double *scratch;
    
	
	scratch = FFTW_MALLOC(local_size_u); // index_u
	
	for (i=0; i<=local_size_u-1; ++i){
		scratch[i] = 0.0;
		Nu_copy[i] = 0.0;
	}
	
	
	// *****************************************************************
	// *****************************************************************
	uu1 = FFTW_MALLOC(local_size_ct); // index_ct
	uu3 = FFTW_MALLOC(local_size_ct); // index_ct
    
	for (i=0; i<=local_size_ct-1; ++i){
		uu1[i] = (9.0*inter->u_phys_ct1[i] - 1.0*inter->u_phys_ct3[i])/8.0\
        *inter->u_phys_ct1[i];
		uu3[i] = (9.0*inter->u_phys_ct1[i] - 1.0*inter->u_phys_ct3[i])/8.0\
        *inter->u_phys_ct3[i];
	}
	
	derivative_x1_ct2u(uu1, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += 9.0/8.0*scratch[i];
	}
    
	// *****************************************************************
	// for LES wall ODE
	const int j_log = par->j_log;
	if (par->les == 1 && par->bc == 1){
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
				index1 = par->index_u[i][j_log][k];
				index2 = par->index_u[i][ny-1-j_log][k];
				index3 = par->index_u[i][j_log+1][k];
				index4 = par->index_u[i][ny-2-j_log][k];
				
				if (par->ode_choice == 1) {
					inter->bc_duudx_b[index] = 9.0/8.0*scratch[index1];
					inter->bc_duudx_t[index] = 9.0/8.0*scratch[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_duudx_b[index]
					= 9.0/8.0*(scratch[index1]+scratch[index3])/2.0;
					inter->bc_duudx_t[index]
					= 9.0/8.0*(scratch[index2]+scratch[index4])/2.0;
				}
			}
		}
	}
	// *****************************************************************
	
	derivative_x3_ct2u(uu3, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += - 1.0/8.0*scratch[i];
	}
	
	// *****************************************************************
	// for LES wall ODE
	if (par->les == 1 && par->bc == 1){
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
				index1 = par->index_u[i][j_log][k];
				index2 = par->index_u[i][ny-1-j_log][k];
				index3 = par->index_u[i][j_log+1][k];
				index4 = par->index_u[i][ny-2-j_log][k];
				
				if (par->ode_choice == 1) {
					inter->bc_duudx_b[index] += - 1.0/8.0*scratch[index1];
					inter->bc_duudx_t[index] += - 1.0/8.0*scratch[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_duudx_b[index]
					+= - 1.0/8.0*(scratch[index1]+scratch[index3])/2.0;
					inter->bc_duudx_t[index]
					+= - 1.0/8.0*(scratch[index2]+scratch[index4])/2.0;
				}
				
			}
		}
        
        
		// transfer into fourier space
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_duudx_b, inter->bc_duudx_b);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_duudx_t, inter->bc_duudx_t);
		for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
			inter->bc_duudx_t[i] = inter->bc_duudx_t[i]/(double)(nz);
			inter->bc_duudx_b[i] = inter->bc_duudx_b[i]/(double)(nz);
		}
	}
	// *****************************************************************
    
	fftw_free(uu1);
	fftw_free(uu3);
	
	
	// *****************************************************************
	// *****************************************************************
	uv1 = FFTW_MALLOC(local_size_cn); // index_cn
	uv3 = FFTW_MALLOC(local_size_cn); // index_cn
	
	for (i=0; i<=local_size_cn-1; ++i){
		uv1[i] = (9.0*inter->v_phys_cn1[i] - 1.0*inter->v_phys_cn3[i])/8.0 \
        *inter->u_phys_cn1[i];
		uv3[i] = (9.0*inter->v_phys_cn1[i] - 1.0*inter->v_phys_cn3[i])/8.0 \
        *inter->u_phys_cn3[i];
	}
	
	// *****************************************
	// additional values of uv3 at j = -1 & ny+1
	index = 0;
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			uv3[par->index_cn[i][-1][k]] \
            = 27.0*uv1[par->index_cn[i][0][k]] \
            - uv3[par->index_cn[i][1][k]] \
            - uv3[par->index_cn[i][0][k]] + 0.0;
            /*+ 0.0 is only for the channel flow.*/
            
			
			uv3[par->index_cn[i][ny+1][k]] \
            = 27.0*uv1[par->index_cn[i][ny][k]] \
            - uv3[par->index_cn[i][ny-1][k]] \
            - uv3[par->index_cn[i][ny][k]] + 0.0;
            /*+ 0.0 is only for the channel flow.*/
            
			
			inter->uv3_bot[index] = uv3[par->index_cn[i][-1][k]];
			inter->uv3_top[index] = uv3[par->index_cn[i][ny+1][k]];
			index += 1;
		}
	}
	// additional values of uv3 at j = -1 & ny+1
	// *****************************************
	
	
	
	derivative_y1_cn2u(uv1, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += 9.0/8.0*scratch[i];
	}
	
	derivative_y3_cn2u(uv3, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(uv1);
	fftw_free(uv3);
	
	// *****************************************************************
	fftw_execute_r2r(planptr->p1d_z_u, Nu_copy, Nu_copy);
    
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] = Nu_copy[i]/(double)nz; // scale
	}
	
	// *****************************************************************
	// *****************************************************************
	uw  = FFTW_MALLOC(local_size_ued); // index_ued
	
	for (i=0; i<=local_size_ued-1; ++i){
		uw[i] = (9.0*inter->w_phys_ued1[i] - 1.0*inter->w_phys_ued3[i])/8.0\
        *inter->u_phys_ued[i];
        
	}
	
	// *****************************************************************
	fftw_execute_r2r(planptr->p1d_z_ued, uw, uw);
    
	for (i=0; i<=local_size_ued-1; ++i){
		uw[i] = uw[i]/(double)nz; // scale
	}
	
	derivative_z_ued2u(uw, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += scratch[i];
	}
	
	// *****************************************************************
	// for LES wall ODE
	if (par->les == 1 && par->bc == 1){
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
				index1 = par->index_u[i][j_log][k];
				index2 = par->index_u[i][ny-1-j_log][k];
				index3 = par->index_u[i][j_log+1][k];
				index4 = par->index_u[i][ny-2-j_log][k];
				
				if (par->ode_choice == 1) {
					inter->bc_duwdz_b[index] = scratch[index1];
					inter->bc_duwdz_t[index] = scratch[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_duwdz_b[index] = (scratch[index1]+scratch[index3])/2.0;
					inter->bc_duwdz_t[index] = (scratch[index2]+scratch[index4])/2.0;
				}
			}
		}
	}
	
	// values are already in fourier space
	// *****************************************************************
	
	
	
	
	fftw_free(uw);
    
	fftw_free(scratch);
    
	return;
}




/* ******************************************************************* */
/* convective (for rhs): given v get Nv(v)							   */
/* v = vector(0, nq_y-1)											   */
/* Nv = vector(0, nq_y-1)											   */
/* ******************************************************************* */
void get_convective_div__4_Nv(
                              interpolated *inter,
                              double *Nv_copy,
                              fftwplans *planptr,
                              parameters *par)
{
    
	int i;
	const int nz = par->nz;
	
	const int local_size_v = par->local_size_v;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	
	double *uv1, *uv3;
	double *vv1, *vv3;
	double *vw;
	double *scratch;
    
	
	scratch = FFTW_MALLOC(local_size_v); // index_wp
	
	for (i=0; i<=local_size_v-1; ++i){
		scratch[i] = 0.0;
		Nv_copy[i] = 0.0;
	}
	
	
	// *****************************************************************
	// *****************************************************************
	uv1 = FFTW_MALLOC(local_size_cn); // index_cn
	uv3 = FFTW_MALLOC(local_size_cn); // index_cn
    
	for (i=0; i<=local_size_cn-1; ++i){
		uv1[i] = (9.0*inter->u_phys_cn1[i] - 1.0*inter->u_phys_cn3[i])/8.0\
        *inter->v_phys_cn1[i];
		uv3[i] = (9.0*inter->u_phys_cn1[i] - 1.0*inter->u_phys_cn3[i])/8.0\
        *inter->v_phys_cn3[i];
	}
	
	
	derivative_x1_cn2v(uv1, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += 9.0/8.0*scratch[i];
	}
    
	derivative_x3_cn2v(uv3, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(uv1);
	fftw_free(uv3);
	
	
	// *****************************************************************
	// *****************************************************************
	vv1 = FFTW_MALLOC(local_size_ct); // index_ct
	vv3 = FFTW_MALLOC(local_size_ct); // index_ct
	
	for (i=0; i<=local_size_ct-1; ++i){
		vv1[i] = (9.0*inter->v_phys_ct1[i] - 1.0*inter->v_phys_ct3[i])/8.0\
        *inter->v_phys_ct1[i];
		vv3[i] = (9.0*inter->v_phys_ct1[i] - 1.0*inter->v_phys_ct3[i])/8.0\
        *inter->v_phys_ct3[i];
	}
	
	derivative_y1_ct2v(vv1, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += 9.0/8.0*scratch[i];
	}
	
	derivative_y3_ct2v(vv3, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(vv1);
	fftw_free(vv3);
	
	// *****************************************************************
	fftw_execute_r2r(planptr->p1d_z_v, Nv_copy, Nv_copy);
    
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] = Nv_copy[i]/(double)nz; // scale
	}
	
	// *****************************************************************
	// *****************************************************************
	vw  = FFTW_MALLOC(local_size_ved); // index_ved
	
	for (i=0; i<=local_size_ved-1; ++i){
		vw[i] = (9.0*inter->w_phys_ved1[i] - 1.0*inter->w_phys_ved3[i])/8.0 \
        *inter->v_phys_ved[i];
        
	}
	
	// *****************************************************************
	fftw_execute_r2r(planptr->p1d_z_ved, vw, vw);
	for (i=0; i<=local_size_ved-1; ++i){
		vw[i] = vw[i]/(double)nz;
	}
	
	derivative_z_ved2v(vw, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += scratch[i];
	}
	
	fftw_free(vw);
	
	fftw_free(scratch);
    
	return;
}




/* ******************************************************************* */
/* convective (for rhs): given w get Nw(w)							   */
/* w = vector(0, nq_w-1)											   */
/* Nw = vector(0, nq_w-1)											   */
/* ******************************************************************* */
void get_convective_div__4_Nw(
                              interpolated *inter,
                              double *Nw_copy,
                              fftwplans *planptr,
                              parameters *par)
{
    
	int i, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int local_size_wp = par->local_size_wp;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	
	
	double *uw1, *uw3;
	double *vw1, *vw3;
	double *ww;
	double *scratch;
    
	
	scratch = FFTW_MALLOC(local_size_wp); // index_wp
	
	for (i=0; i<=local_size_wp-1; ++i){
		scratch[i] = 0.0;
		Nw_copy[i] = 0.0;
	}
	
	
	// *****************************************************************
	// *****************************************************************
	uw1 = FFTW_MALLOC(local_size_ued); // index_ued
	uw3 = FFTW_MALLOC(local_size_ued); // index_ued
    
	for (i=0; i<=local_size_ued-1; ++i){
		uw1[i] = inter->u_phys_ued[i]*inter->w_phys_ued1[i];
		uw3[i] = inter->u_phys_ued[i]*inter->w_phys_ued3[i];
	}
	
	
	derivative_x1_ued2wp(uw1, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += 9.0/8.0*scratch[i];
	}
    
	derivative_x3_ued2wp(uw3, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(uw1);
	fftw_free(uw3);
	
	
	// *****************************************************************
	// *****************************************************************
	vw1 = FFTW_MALLOC(local_size_ved); // index_ved
	vw3 = FFTW_MALLOC(local_size_ved); // index_ved
	
	for (i=0; i<=local_size_ved-1; ++i){
		vw1[i] = inter->v_phys_ved[i]*inter->w_phys_ved1[i];
		vw3[i] = inter->v_phys_ved[i]*inter->w_phys_ved3[i];
	}
	
	
	// *****************************************
	// additional values of vw3 at j = -1 & ny+1
	index = 0;
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			vw3[par->index_ved[i][-1][k]] \
            = 27.0*vw1[par->index_ved[i][0][k]] \
            - vw3[par->index_ved[i][1][k]] \
            - vw3[par->index_ved[i][0][k]] + 0.0;
            /*+ 0.0 is only for the channel flow.*/
            
			
			vw3[par->index_ved[i][ny+1][k]] \
            = 27.0*vw1[par->index_ved[i][ny][k]] \
            - vw3[par->index_ved[i][ny-1][k]] \
            - vw3[par->index_ved[i][ny][k]] + 0.0;
            /*+ 0.0 is only for the channel flow.*/
            
			
			inter->wv3_bot[index] = vw3[par->index_ved[i][-1][k]];
			inter->wv3_top[index] = vw3[par->index_ved[i][ny+1][k]];
			index += 1;
		}
	}
	// additional values of vw3 at j = -1 & ny+1
	// *****************************************
	
	derivative_y1_ved2wp(vw1, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += 9.0/8.0*scratch[i];
	}
	
	derivative_y3_ved2wp(vw3, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(vw1);
	fftw_free(vw3);
	
	// *****************************************************************
	fftw_execute_r2r(planptr->p1d_z_wp, Nw_copy, Nw_copy);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] = Nw_copy[i]/(double)nz; // scale
	}
	
	// *****************************************************************
	// *****************************************************************
	ww  = FFTW_MALLOC(local_size_ct); // index_ct
	
	for (i=0; i<=local_size_ct-1; ++i){
		ww[i] = inter->w_phys_ct[i]*inter->w_phys_ct[i];
	}
	
	// *****************************************************************
	fftw_execute_r2r(planptr->p1d_z_ct, ww, ww);
	for (i=0; i<=local_size_ct-1; ++i){
		ww[i] = ww[i]/(double)nz; // scale
	}
	
	derivative_z_ct2wp(ww, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += scratch[i];
	}
	
	fftw_free(ww);
	
	
	fftw_free(scratch);
    
	return;
}
