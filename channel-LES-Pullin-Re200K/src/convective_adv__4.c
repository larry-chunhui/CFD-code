/*
 *  convective_adv__4.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Calculates the advective form of the convective terms with
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
#include "differentiate.h"
#include "interpolation__4.h"
#include "convective_div__4.h"
#include "convective_adv__4.h"


#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

/* ******************************************************************* */
/* convective adi form (for rhs): given u get Nu(u)					   */
/* u = vector(0, nq_x-1)											   */
/* Nu = vector(0, nq_x-1)											   */
/* v = vector(0, nq_y-1)											   */
/* Nv = vector(0, nq_y-1)											   */
/* w = vector(0, nq_w-1)											   */
/* Nw = vector(0, nq_w-1)											   */
/* ******************************************************************* */
void get_convective_adv__4(
                           interpolated *inter,
                           double *Nu_copy,
                           double *Nv_copy,
                           double *Nw_copy,
                           fftwplans *planptr,
                           parameters *par)
{
	
	get_convective_adv__4_Nu(inter, Nu_copy, planptr, par);
	get_convective_adv__4_Nv(inter, Nv_copy, planptr, par);
	get_convective_adv__4_Nw(inter, Nw_copy, planptr, par);
    
	return;
}

/* ******************************************************************* */
/* convective (for rhs): given u get Nu(u)							   */
/* u = vector(0, nq_x-1)											   */
/* Nu = vector(0, nq_x-1)											   */
/* ******************************************************************* */
void get_convective_adv__4_Nu(
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
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	const double dy = par->dy;
	
	double *udu1, *udu3;
	double *vdu1, *vdu3;
	double *w_x3, *dudz;
	double *scratch;
    
	
	scratch = FFTW_MALLOC(local_size_u); // index_u
	
	for (i=0; i<=local_size_u-1; ++i){
		scratch[i] = 0.0;
		Nu_copy[i] = 0.0;
	}
	
	
	udu1 = FFTW_MALLOC(local_size_ct); // index_ct
	udu3 = FFTW_MALLOC(local_size_ct); // index_ct
	
	// dudx1
	derivative_x1_ued2ct__4(inter->u_phys_ued, udu1, par);
	// udu1
	for (i=0; i<=local_size_ct-1; ++i){
		udu1[i] = (9.0*inter->u_phys_ct1[i] - 1.0*inter->u_phys_ct3[i])/8.0\
        *udu1[i];
	}
	
	// dudx3
	derivative_x3_ued2ct__4(inter->u_phys_ued, udu3, par);
	// udu3
	for (i=0; i<=local_size_ct-1; ++i){
		udu3[i] = (9.0*inter->u_phys_ct1[i] - 1.0*inter->u_phys_ct3[i])/8.0\
        *udu3[i];
	}
	
	
	interpolate_x1_ct2u__4(udu1, scratch, par);
	
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
					inter->bc_ududx_b[index] = 9.0/8.0*scratch[index1];
					inter->bc_ududx_t[index] = 9.0/8.0*scratch[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_ududx_b[index]
					= 9.0/8.0*(scratch[index1]+scratch[index3])/2.0;
					inter->bc_ududx_t[index]
					= 9.0/8.0*(scratch[index2]+scratch[index4])/2.0;
				}
			}
		}
	}
	// *****************************************************************
	
	
	
	interpolate_x3_ct2u__4(udu3, scratch, par);
	
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
					inter->bc_ududx_b[index] += - 1.0/8.0*scratch[index1];
					inter->bc_ududx_t[index] += - 1.0/8.0*scratch[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_ududx_b[index]
					+= - 1.0/8.0*(scratch[index1]+scratch[index3])/2.0;
					inter->bc_ududx_t[index]
					+= - 1.0/8.0*(scratch[index2]+scratch[index4])/2.0;
				}
			}
		}
        fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_ududx_b, inter->bc_ududx_b);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_ududx_t, inter->bc_ududx_t);
		for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
			inter->bc_ududx_t[i] = inter->bc_ududx_t[i]/(double)(nz);
			inter->bc_ududx_b[i] = inter->bc_ududx_b[i]/(double)(nz);
		}
	}
	// *****************************************************************
	
	
	fftw_free(udu1);
	fftw_free(udu3);
	
	vdu1 = FFTW_MALLOC(local_size_cn); // index_cn
	vdu3 = FFTW_MALLOC(local_size_cn); // index_cn
	
	// dudy1
	derivative_y1_ued2cn__4(inter->u_phys_ued, vdu1, par);
	// vdu1
	for (i=0; i<=local_size_cn-1; ++i){
		vdu1[i] = (9.0*inter->v_phys_cn1[i] - 1.0*inter->v_phys_cn3[i])/8.0\
        *vdu1[i];
	}
	
	// dudy3
	derivative_y3_ued2cn__4(inter->u_phys_ued, vdu3, par);
    
	// vdu3
	for (i=0; i<=local_size_cn-1; ++i){
		vdu3[i] = (9.0*inter->v_phys_cn1[i] - 1.0*inter->v_phys_cn3[i])/8.0\
        *vdu3[i];
	}
	
	// additional values of uv3 at j = -1 & ny+1
	index = 0;
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			vdu3[par->index_cn[i][-1][k]] \
            = ( 2.0*inter->u_phys_ued[par->index_ued[i][0][k]] \
               *(9.0*inter->v_phys_cn1[par->index_cn[i][-1][k]] - 1.0*inter->v_phys_cn3[par->index_cn[i][-1][k]])/8.0 \
               - 2.0*inter->uv3_bot[index] )/(3.0*dy);
            
            
			
			vdu3[par->index_cn[i][ny+1][k]] \
            = ( - 2.0*inter->u_phys_ued[par->index_ued[i][ny-1][k]] \
               *(9.0*inter->v_phys_cn1[par->index_cn[i][ny+1][k]] - 1.0*inter->v_phys_cn3[par->index_cn[i][ny+1][k]])/8.0 \
               + 2.0*inter->uv3_top[index] )/(3.0*dy);
            
            
			index += 1;
		}
	}
	
	interpolate_y1_cn2u__4(vdu1, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += 9.0/8.0*scratch[i];
	}
	
	interpolate_y3_cn2u__4(vdu3, scratch, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(vdu1);
	fftw_free(vdu3);
	
    
	
	w_x3 = FFTW_MALLOC(local_size_u); // index_u
	dudz  = FFTW_MALLOC(local_size_u); // index_u
	
	// dudz
	derivative_z_ued2u(inter->u_ued, dudz, par);
	fftw_execute_r2r(planptr->p1d_invz_u, dudz, dudz);
	// w_x1
	interpolate_x1_ct2u__4(inter->w_phys_ct, scratch, par);
	// w_x3
	interpolate_x3_ct2u__4(inter->w_phys_ct, w_x3, par);
	// calculate wdu and add
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += (9.0*scratch[i] - 1.0*w_x3[i])/8.0 \
        *dudz[i];
		scratch[i]  = (9.0*scratch[i] - 1.0*w_x3[i])/8.0 \
        *dudz[i];
	}
	
	
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
					inter->bc_wdudz_b[index] = scratch[index1];
					inter->bc_wdudz_t[index] = scratch[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_wdudz_b[index]
					= (scratch[index1]+scratch[index3])/2.0;
					inter->bc_wdudz_t[index]
					= (scratch[index2]+scratch[index4])/2.0;
				}
			}
		}
        
        fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_wdudz_b, inter->bc_wdudz_b);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_wdudz_t, inter->bc_wdudz_t);
		for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
			inter->bc_wdudz_t[i] = inter->bc_wdudz_t[i]/(double)(nz);
			inter->bc_wdudz_b[i] = inter->bc_wdudz_b[i]/(double)(nz);
		}
	}
	
	
	fftw_free(w_x3);
	fftw_free(dudz);
	fftw_free(scratch);
    
	fftw_execute_r2r(planptr->p1d_z_u, Nu_copy, Nu_copy);
    
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] = Nu_copy[i]/(double)nz; // scale
	}
    
	return;
}



/* ******************************************************************* */
/* convective (for rhs): given v get Nv(v)							   */
/* v = vector(0, nq_y-1)											   */
/* Nv = vector(0, nq_y-1)											   */
/* ******************************************************************* */
void get_convective_adv__4_Nv(
                              interpolated *inter,
                              double *Nv_copy,
                              fftwplans *planptr,
                              parameters *par)
{
	
	int i;
	const int nz = par->nz;
	const int local_size_v = par->local_size_v;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	
	double *udv1, *udv3;
	double *vdv1, *vdv3;
	double *w_y3, *dvdz;
	double *scratch;
    
	
	scratch = FFTW_MALLOC(local_size_v); // index_u
	
	for (i=0; i<=local_size_v-1; ++i){
		scratch[i] = 0.0;
		Nv_copy[i] = 0.0;
	}
	
    udv1 = FFTW_MALLOC(local_size_cn); // index_cn
	udv3 = FFTW_MALLOC(local_size_cn); // index_cn
	
	// dvdx1
	derivative_x1_ved2cn__4(inter->v_phys_ved, udv1, par);
	// udv1
	for (i=0; i<=local_size_cn-1; ++i){
		udv1[i] = (9.0*inter->u_phys_cn1[i] - 1.0*inter->u_phys_cn3[i])/8.0\
        *udv1[i];
	}
	
	// get dvdx3
	derivative_x3_ved2cn__4(inter->v_phys_ved, udv3, par);
	// get udu3
	for (i=0; i<=local_size_cn-1; ++i){
		udv3[i] = (9.0*inter->u_phys_cn1[i] - 1.0*inter->u_phys_cn3[i])/8.0\
        *udv3[i];
	}
	
	
	interpolate_x1_cn2v__4(udv1, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += 9.0/8.0*scratch[i];
	}
    
	interpolate_x3_cn2v__4(udv3, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(udv1);
	fftw_free(udv3);
	
	vdv1 = FFTW_MALLOC(local_size_ct); // index_ct
	vdv3 = FFTW_MALLOC(local_size_ct); // index_ct
	
	// dvdy1
	derivative_y1_ved2ct__4(inter->v_phys_ved, vdv1, par);
	// vdv1
	for (i=0; i<=local_size_ct-1; ++i){
		vdv1[i] = (9.0*inter->v_phys_ct1[i] - 1.0*inter->v_phys_ct3[i])/8.0\
        *vdv1[i];
	}
	
	// dvdy3
	derivative_y3_ved2ct__4(inter->v_phys_ved, vdv3, par);
	// vdv3
	for (i=0; i<=local_size_ct-1; ++i){
		vdv3[i] = (9.0*inter->v_phys_ct1[i] - 1.0*inter->v_phys_ct3[i])/8.0\
        *vdv3[i];
	}
    
	
	interpolate_y1_ct2v__4(vdv1, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += 9.0/8.0*scratch[i];
	}
	
	interpolate_y3_ct2v__4(vdv3, scratch, par);
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(vdv1);
	fftw_free(vdv3);
	
    
	
	w_y3  = FFTW_MALLOC(local_size_v); // index_u
	dvdz  = FFTW_MALLOC(local_size_v); // index_u
	
	// dvdz
	derivative_z_ved2v(inter->v_ved, dvdz, par);
	fftw_execute_r2r(planptr->p1d_invz_v, dvdz, dvdz);
	// w_y1
	interpolate_y1_ct2v__4(inter->w_phys_ct, scratch, par);
	// w_y3
	interpolate_y3_ct2v__4(inter->w_phys_ct, w_y3, par);
	// calculate wdu and add
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += (9.0*scratch[i] - 1.0*w_y3[i])/8.0 \
        *dvdz[i];
	}
	
	fftw_free(w_y3);
	fftw_free(dvdz);
    
	fftw_free(scratch);
    
	fftw_execute_r2r(planptr->p1d_z_v, Nv_copy, Nv_copy);
    
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] = Nv_copy[i]/(double)nz; // scale
	}
    
	return;
}


/* ******************************************************************* */
/* convective (for rhs): given w get Nw(w)							   */
/* v = vector(0, nq_z-1)											   */
/* Nv = vector(0, nq_z-1)											   */
/* ******************************************************************* */
void get_convective_adv__4_Nw(
                              interpolated *inter,
                              double *Nw_copy,
                              fftwplans *planptr,
                              parameters *par)
{
	
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int local_size_wp = par->local_size_wp;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const double dy = par->dy;
	
	double *udw1, *udw3;
	double *vdw1, *vdw3;
	double *wdw, *dwdz;
	double *scratch;
    
	
	scratch = FFTW_MALLOC(local_size_wp); // index_u
	
	for (i=0; i<=local_size_wp-1; ++i){
		scratch[i] = 0.0;
		Nw_copy[i] = 0.0;
	}
	
	
	udw1 = FFTW_MALLOC(local_size_ued); // index_ued
	udw3 = FFTW_MALLOC(local_size_ued); // index_ued
	
	// dwdx1
	derivative_x1_ct2ued__4(inter->w_phys_ct, udw1, par);
	// udw1
	for (i=0; i<=local_size_ued-1; ++i){
		udw1[i] *= inter->u_phys_ued[i];
        
	}
	
	// dwdx3
	derivative_x3_ct2ued__4(inter->w_phys_ct, udw3, par);
	// udw3
	for (i=0; i<=local_size_ued-1; ++i){
		udw3[i] *= inter->u_phys_ued[i];
        
	}
	
	
	interpolate_x1_ued2wp__4(udw1, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += 9.0/8.0*scratch[i];
	}
    
	interpolate_x3_ued2wp__4(udw3, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(udw1);
	fftw_free(udw3);
	
	vdw1 = FFTW_MALLOC(local_size_ved); // index_ved
	vdw3 = FFTW_MALLOC(local_size_ved); // index_ved
	
	// dwdy1
	derivative_y1_ct2ved__4(inter->w_phys_ct, vdw1, par);
	// vdw1
	for (i=0; i<=local_size_ved-1; ++i){
		vdw1[i] *= inter->v_phys_ved[i];
        
	}
	
	// dwdy3
	derivative_y3_ct2ved__4(inter->w_phys_ct, vdw3, par);
    
	// vdw3
	for (i=0; i<=local_size_ved-1; ++i){
		vdw3[i] *= inter->v_phys_ved[i];
        
	}
	
	
	
	// additional values of uv3 at j = -1 & ny+1
	index = 0;
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			vdw3[par->index_ved[i][-1][k]] \
            = ( 2.0*inter->w_phys_ct[par->index_ct[i][0][k]]*inter->v_phys_ved[par->index_ved[i][-1][k]] \
               - 2.0*inter->wv3_bot[index] )/(3.0*dy);
            
            
			
			vdw3[par->index_ved[i][ny+1][k]] \
            = ( - 2.0*inter->w_phys_ct[par->index_ct[i][ny-1][k]]*inter->v_phys_ved[par->index_ved[i][ny+1][k]] \
               + 2.0*inter->wv3_top[index] )/(3.0*dy);
            
            
			index += 1;
		}
	}
	
	
	
	
	
	interpolate_y1_ved2wp__4(vdw1, scratch, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += 9.0/8.0*scratch[i];
	}
    
	interpolate_y3_ved2wp__4(vdw3, scratch, par);
    
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += - 1.0/8.0*scratch[i];
	}
    
	fftw_free(vdw1);
	fftw_free(vdw3);
	
    wdw  = FFTW_MALLOC(local_size_wp); // index_wp
	dwdz  = FFTW_MALLOC(local_size_wp); // index_wp
	
	// dwdz
	derivative_z_ct2wp(inter->w_ct, dwdz, par);
	fftw_execute_r2r(planptr->p1d_invz_wp, dwdz, dwdz);
	// calculate wdw and add
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Nw_copy[par->index_wp[i][j][k]] \
                += inter->w_phys_ct[par->index_ct[i][j][k]]*dwdz[par->index_wp[i][j][k]];
			}
		}
	}
	
	
	fftw_free(wdw);
	fftw_free(dwdz);
    
	fftw_free(scratch);
    
	fftw_execute_r2r(planptr->p1d_z_wp, Nw_copy, Nw_copy);
    
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] = Nw_copy[i]/(double)nz; // scale
	}
    
	return;
}



