/*
 *  convective_div__2.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Calculates the divergence form of the convective terms with
 * second order schemes.
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
#include "convective_div__2.h"
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
void get_convective_div__2(
                           interpolated *inter,
                           double *Nu_copy,
                           double *Nv_copy,
                           double *Nw_copy,
                           fftwplans *planptr,
                           parameters *par)
{
	int i, k, index;
	int index1, index2, index3, index4;
	const int ny = par->ny;
	const int nz = par->nz;
	
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const int local_size_u = par->local_size_u;
	const int local_size_v = par->local_size_v;
	const int local_size_wp = par->local_size_wp;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	double *uu, *vv, *ww;
	double *uv, *uw, *vw;
	double *scratch_u, *scratch_v, *scratch_wp;
    
	uu = FFTW_MALLOC(local_size_ct); // index_ct
	vv = FFTW_MALLOC(local_size_ct); // index_ct
	ww = FFTW_MALLOC(local_size_ct); // index_ct
	uv = FFTW_MALLOC(local_size_cn); // index_cn
	uw = FFTW_MALLOC(local_size_ued); // index_ued
	vw = FFTW_MALLOC(local_size_ved); // index_ved
    
	scratch_u = FFTW_MALLOC(local_size_u); // index_u
	scratch_v = FFTW_MALLOC(local_size_v); // index_v
	scratch_wp = FFTW_MALLOC(local_size_wp); // index_wp
	
	for (i=0; i<=local_size_ct-1; ++i){
		uu[i] = 0.0;
		vv[i] = 0.0;
		ww[i] = 0.0;
	}
	
	for (i=0; i<=local_size_cn-1; ++i){
		uv[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ued-1; ++i){
		uw[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ved-1; ++i){
		vw[i] = 0.0;
	}
	
    
	
	
	/* (for example when my_rank == 0)
     uu[], vv[], ww[] is defined at cell center  [start-1->end][-1->ny][0->nz-1]
     index_ct[][][] is used as a location indicator
     
     uv[] is defined at corners of each cell [start->end][0->ny][0->nz-1] (boundary included)
     index_cn[][][] is used as a location indicator
     
     uw[] is defined at edges of each cell [start-1->end][-1->ny][0->nz-1] (boundary included)
     similar to location of velocity u[]; [start->local_nx-1][0->ny-1][0->nz-1]
     index_ued[][][] is used as a location indicator
     
     vw[] is defined at edges of each cell [start-1->end][0->ny][0->nz-1] (boundary included)
     similar to location of velocity v[]; [start->end-1][1->ny-1][0->nz-1]
     index_ved[][][] is used as a location indicator
     */
	
	form_product_div__2(inter, uu, vv, ww, uv, uw, vw, par, planptr);
    
	
	/***********************************************************************************
     d(uw)/dz defined at the same position is velocity u[] (edge)
     index_u[][][] is used as a location indicator
     while index_ued[][][] is for uw[]
     */
	
	/***********************************************************************************
     d(vw)/dz defined at the same position is velocity v[] (edge)
     index_v[][][] is used as a location indicator
     while index_ved[][][] is for vw[]
     */
	
	/***********************************************************************************
     d(ww)/dz defined at the same position is velocity w[] (center)
     index_wp[][][] is used as a location indicator
     same for ww[]
     */
	
	//***********************************************************************************
	/***********************************************************************************
     Nu = d(uu)/dx + d(uv)/dy + d(uw)/dz defined at the same as velocity u[]
     uu -> index_ct[][][]
     uv -> index_cn[][][]
     d(uw)/dz -> index_u[][][]
     */
	/***********************************************************************************
     Nv = d(uv)/dx + d(vv)/dy + d(vw)/dz, defined at the same as velocity v[]
     uv -> index_cn[][][]
     vv -> index_ct[][][]
     d(vw)/dz -> index_v[][][]
     */
	
	/************************************************************************************
     Nw = d(uw)/dx + d(vw)/dy + d(ww)/dz, defined at the same as velocity w[]
     uw -> index_ued[][][]
     vw -> index_ved[][][]
     d(ww)/dz -> index_wp[][][]
     */
    
    
	
	
	/***********************************************************************************/
	/***********************************************************************************/
	derivative_x1_ct2u(uu, scratch_u, par);
	derivative_x1_cn2v(uv, scratch_v, par);
	derivative_x1_ued2wp(uw, scratch_wp, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] = scratch_u[i];
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
					inter->bc_duudx_b[index] = scratch_u[index1];
					inter->bc_duudx_t[index] = scratch_u[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_duudx_b[index]
					= (scratch_u[index1]+scratch_u[index3])/2.0;
					inter->bc_duudx_t[index]
					= (scratch_u[index2]+scratch_u[index4])/2.0;
				}
			}
		}
	}
	// *****************************************************************
	
	
	
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] = scratch_v[i];
	}
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] = scratch_wp[i];
	}
	
	
	/***********************************************************************************/
	/***********************************************************************************/
	derivative_y1_cn2u(uv, scratch_u, par);
	derivative_y1_ct2v(vv, scratch_v, par);
	derivative_y1_ved2wp(vw, scratch_wp, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += scratch_u[i];
	}
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += scratch_v[i];
	}
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += scratch_wp[i];
	}
	
	
	/***********************************************************************************/
	/***********************************************************************************/
	derivative_z_ued2u(uw, scratch_u, par);
	derivative_z_ved2v(vw, scratch_v, par);
	derivative_z_ct2wp(ww, scratch_wp, par);
	
	for (i=0; i<=local_size_u-1; ++i){
		Nu_copy[i] += scratch_u[i];
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
					inter->bc_duwdz_b[index] = scratch_u[index1];
					inter->bc_duwdz_t[index] = scratch_u[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_duwdz_b[index] = (scratch_u[index1]+scratch_u[index3])/2.0;
					inter->bc_duwdz_t[index] = (scratch_u[index2]+scratch_u[index4])/2.0;
				}
			}
		}
	}
	
	// values are already in fourier space
	// *****************************************************************
	
	
	for (i=0; i<=local_size_v-1; ++i){
		Nv_copy[i] += scratch_v[i];
	}
	
	for (i=0; i<=local_size_wp-1; ++i){
		Nw_copy[i] += scratch_wp[i];
	}
	
	
    
	// free allocated memory
	fftw_free(uu);
	fftw_free(vv);
	fftw_free(ww);
	fftw_free(uv);
	fftw_free(uw);
	fftw_free(vw);
	
	fftw_free(scratch_u);
	fftw_free(scratch_v);
	fftw_free(scratch_wp);
	
	
	return;
}



/* ******************************************************************* */
/* Basic procedure
 
 interpolate u at center (u_ct) and corner (u_cn)
 interpolate v at center (v_ct) and corner (v_cn)
 interpolate w at edge (w_ed)
 
 u -> u_ed, u_ct, u_cn
 v -> v_ed, v_ct, v_cn
 w -> w, w_ued, w_ved
 
 inverse fft in z-direction
 multiply to form
 Nu; u_ct*u_ct, u*w_ed, u_cn*v_cn
 Nv; v_ct*v_ct, v*w_ed, u_cn*v_cn
 Nw; w*w, u*w_ued, v*w_ved
 */
/* ******************************************************************* */
void form_product_div__2(
                         interpolated *inter,
                         double *uu,
                         double *vv,
                         double *ww,
                         double *uv,
                         double *uw,
                         double *vw,
                         parameters *par,
                         fftwplans *planptr)
{
	
	int i;
	const int nz = par->nz;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
    
	// calculate necessary products at each location
	calc_multiply_ct(inter->u_phys_ct1, inter->u_phys_ct1, uu, par); // uu defined at center
	calc_multiply_ct(inter->v_phys_ct1, inter->v_phys_ct1, vv, par); // vv defined at center
	calc_multiply_ct(inter->w_phys_ct, inter->w_phys_ct, ww, par); // ww defined at center
	calc_multiply_ued(inter->u_phys_ued, inter->w_phys_ued1, uw, par); // uw defuned at edge
	calc_multiply_ved(inter->v_phys_ved, inter->w_phys_ved1, vw, par); // vw defined at edge
	calc_multiply_cn(inter->u_phys_cn1, inter->v_phys_cn1, uv, par); // uv defined at corner
	
	
	// forward fft back into wave space
	fftw_execute_r2r(planptr->p1d_z_ct, uu, uu);
	fftw_execute_r2r(planptr->p1d_z_ct, vv, vv);
	fftw_execute_r2r(planptr->p1d_z_ct, ww, ww);
	// index_ct[][][];[start-1->end][-1->ny][0->nz-1]
	for (i=0; i<=local_size_ct-1; ++i){
		uu[i] = uu[i]/(double)nz;
		vv[i] = vv[i]/(double)nz;
		ww[i] = ww[i]/(double)nz; // scale
	}
	
	fftw_execute_r2r(planptr->p1d_z_cn, uv, uv);
	// index_cn[][][]; [start->end][0->ny][0->nz-1]
	for (i=0; i<=local_size_cn-1; ++i){
		uv[i] = uv[i]/(double)nz;
	}
	
	fftw_execute_r2r(planptr->p1d_z_ued, uw, uw);
	// index_ued[][][]; [start-1->end][-1->ny][0->nz-1]
	for (i=0; i<=local_size_ued-1; ++i){
		uw[i] = uw[i]/(double)nz;
	}
	
	fftw_execute_r2r(planptr->p1d_z_ved, vw, vw);
	// index_ved[][][]; [start-1->end][0->ny][0->nz-1]
	for (i=0; i<=local_size_ved-1; ++i){
		vw[i] = vw[i]/(double)nz;
	}
	
	return;
}

