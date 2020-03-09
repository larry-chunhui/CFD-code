/*
 * velgrad_tensor.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Calculates the velocity gradients at cn.
 * They are used as some of the inputs to the stretched vortex model.
 *
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
#include "interpolation.h"
#include "interpolation__2.h"
#include "interpolation__4.h"
#include "differentiate.h"
#include "velgrad_tensor.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


void get_dwdy_ved_wall_ver2(double *w, double *dwdy, parameters *par);

/* ******************************************************************* */
/* get velocity gradient tensor defined at corner					   */
/* these values will be used to calculate subgrid stress tensor		   */
/* Originally,														   */
/* dudx, dvdy and dwdz are defined at center						   */
/* will be interpolated to v-edge, then to corner					   */
/* dudz, dwdx are defined at u-edge, will be interpolated to corner    */
/* dvdz, dwdy are defined at v-edge, will be interpolated to corner    */
/* dudy, dvdx are defined at corner (no interpolation needed)		   */
/* ******************************************************************* */
void get_velgrad_tensor_cn(
						   variables *var,
						   interpolated *inter,
						   parameters *par,
						   fftwplans *planptr)
{
	
	int i, k;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *scratch1;
	double *scratch2;
	double *scratch3;
	
	double *u_b, *u_t;
	double *w_b, *w_t;
	
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	
	u_b = FFTW_MALLOC((nz)*(local_nx+1));
	u_t = FFTW_MALLOC((nz)*(local_nx+1));
	w_b = FFTW_MALLOC((nz)*(local_nx+1));
	w_t = FFTW_MALLOC((nz)*(local_nx+1));
	
	
	scratch1 = FFTW_MALLOC(local_size_ct);
	scratch2 = FFTW_MALLOC(local_size_ct);
	
	
	// get dudx *******************************************
    // from ued onto ct
	derivative_x1_ued2ct__4(inter->u_phys_ued, scratch1, par);
	derivative_x3_ued2ct__4(inter->u_phys_ued, scratch2, par);
	
	// u_ct is not defined at i = start-3, end+2, end+1
	for (i=0; i<=local_size_ct-1; ++i){
		scratch1[i] = ( 9.0*scratch1[i] - 1.0*scratch2[i] )/8.0;
	}
	// from ct to cn
	interpolate_ct2cn__4(scratch1, var->dudx, par);
	
	
	
	// get dvdy *******************************************
	derivative_y1_ved2ct__4(inter->v_phys_ved, scratch1, par);
	derivative_y3_ved2ct__4(inter->v_phys_ved, scratch2, par);
	
	for (i=0; i<=local_size_ct-1; ++i){
		scratch1[i] = ( 9.0*scratch1[i] - 1.0*scratch2[i] )/8.0;
	}
	
	interpolate_ct2cn__4(scratch1, var->dvdy, par);
	
	// in addition, dvdy at the wall
	{
		
		double *dvdy_ved;
		double *scratch_cn1;
		double *scratch_cn3;
		dvdy_ved = FFTW_MALLOC(local_size_ved);
		scratch_cn1 = FFTW_MALLOC(local_size_cn);
		scratch_cn3 = FFTW_MALLOC(local_size_cn);
		
		for (i=0; i<=local_size_ved-1; ++i){
			dvdy_ved[i] = 0.0;
		}
		
		for (i=start-3; i<=end+2; ++i){
			for (k=0; k<=nz-1; ++k){
				dvdy_ved[par->index_ved[i][0][k]]
				= (-   11.0*inter->v_phys_ved[par->index_ved[i][0][k]]
				   +   18.0*inter->v_phys_ved[par->index_ved[i][1][k]]
				   -    9.0*inter->v_phys_ved[par->index_ved[i][2][k]]
				   +    2.0*inter->v_phys_ved[par->index_ved[i][3][k]])/(6.0*par->dy);
				
				dvdy_ved[par->index_ved[i][ny][k]]
				= (+   11.0*inter->v_phys_ved[par->index_ved[i][ny][k]]
				   -   18.0*inter->v_phys_ved[par->index_ved[i][ny-1][k]]
				   +    9.0*inter->v_phys_ved[par->index_ved[i][ny-2][k]]
				   -    2.0*inter->v_phys_ved[par->index_ved[i][ny-3][k]])/(6.0*par->dy);
				
				dvdy_ved[par->index_ved[i][1][k]]
				= (-    2.0*inter->v_phys_ved[par->index_ved[i][0][k]]
				   -    3.0*inter->v_phys_ved[par->index_ved[i][1][k]]
				   +    6.0*inter->v_phys_ved[par->index_ved[i][2][k]]
				   -    1.0*inter->v_phys_ved[par->index_ved[i][3][k]])/(6.0*par->dy);
				
				dvdy_ved[par->index_ved[i][ny-1][k]]
				= (+    2.0*inter->v_phys_ved[par->index_ved[i][ny][k]]
				   +    3.0*inter->v_phys_ved[par->index_ved[i][ny-1][k]]
				   -    6.0*inter->v_phys_ved[par->index_ved[i][ny-2][k]]
				   +    1.0*inter->v_phys_ved[par->index_ved[i][ny-3][k]])/(6.0*par->dy);
			}
		}
		
		interpolate_x1_ved2cn__4(dvdy_ved, scratch_cn1, par);
		interpolate_x3_ved2cn__4(dvdy_ved, scratch_cn3, par);
		
		for (i=start-1; i<=end+1; ++i){
			for (k=0; k<=nz-1; ++k){
				var->dvdy[par->index_cn[i][0][k]]
				= 9.0/8.0*scratch_cn1[par->index_cn[i][0][k]]
				- 1.0/8.0*scratch_cn3[par->index_cn[i][0][k]];
				
				var->dvdy[par->index_cn[i][1][k]]
				= 9.0/8.0*scratch_cn1[par->index_cn[i][1][k]]
				- 1.0/8.0*scratch_cn3[par->index_cn[i][1][k]];
				
				var->dvdy[par->index_cn[i][ny][k]]
				= 9.0/8.0*scratch_cn1[par->index_cn[i][ny][k]]
				- 1.0/8.0*scratch_cn3[par->index_cn[i][ny][k]];
				
				var->dvdy[par->index_cn[i][ny-1][k]]
				= 9.0/8.0*scratch_cn1[par->index_cn[i][ny-1][k]]
				- 1.0/8.0*scratch_cn3[par->index_cn[i][ny-1][k]];
			}
		}
		
		
		fftw_free(dvdy_ved);
		fftw_free(scratch_cn1);
		fftw_free(scratch_cn3);
		
	}
	
	
	// get dwdz *******************************************
	derivative_z_ct2ct(inter->w_ct, scratch1, par);
	fftw_execute_r2r(planptr->p1d_invz_ct, scratch1, scratch1);
	// from ct to cn
	interpolate_ct2cn__4(scratch1, var->dwdz, par);
	
	fftw_free(scratch1);
	fftw_free(scratch2);
	
	// get dudz ******************************************
    scratch1 = FFTW_MALLOC(local_size_ued);
	scratch2 = FFTW_MALLOC(local_size_ued);
	scratch3 = FFTW_MALLOC(local_size_cn);
    
    // take the derivative of u (ued) w.r.t. z at ued
	derivative_z_ued2ued(inter->u_ued, scratch1, par);
    
    // fft to the physical space
	fftw_execute_r2r(planptr->p1d_invz_ued, scratch1, scratch1);
	
    // interpolate from ued onto cn
	interpolate_y1_ued2cn__4(scratch1, var->dudz, par);
	interpolate_y3_ued2cn__4(scratch1, scratch3, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->dudz[i] = ( 9.0*var->dudz[i] - 1.0*scratch3[i] )/8.0;
	}
	
	
	
	// get dwdx ******************************************
	derivative_x1_ct2ued__4(inter->w_phys_ct, scratch1, par);
	derivative_x3_ct2ued__4(inter->w_phys_ct, scratch2, par);
	
	for (i=0; i<=local_size_ued-1; ++i){
		scratch1[i] = ( 9.0*scratch1[i] - 1.0*scratch2[i] )/8.0;
	}
	
	
	interpolate_y1_ued2cn__4(scratch1, var->dwdx, par);
	interpolate_y3_ued2cn__4(scratch1, scratch3, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->dwdx[i] = ( 9.0*var->dwdx[i] - 1.0*scratch3[i] )/8.0;
	}
    
	fftw_free(scratch1);
	fftw_free(scratch2);
	fftw_free(scratch3);
	
	
	
	
	
	// get dvdz ******************************************
	scratch1 = FFTW_MALLOC(local_size_ved);
	scratch2 = FFTW_MALLOC(local_size_ved);
	scratch3 = FFTW_MALLOC(local_size_cn);
	
    derivative_z_ved2ved(inter->v_ved, scratch1, par);
	
	fftw_execute_r2r(planptr->p1d_invz_ved, scratch1, scratch1);
	
	interpolate_x1_ved2cn__4(scratch1, var->dvdz, par);
	interpolate_x3_ved2cn__4(scratch1, scratch3, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->dvdz[i] = ( 9.0*var->dvdz[i] - 1.0*scratch3[i] )/8.0;
	}
	
	
	// get dwdy ******************************************
	// need some supplimental work at j = 0, ny-1
	// since dwdy at j = 0, ny-1 are set to be zero
	// using the following methods.
	derivative_y1_ct2ved__4(inter->w_phys_ct, scratch1, par);
	derivative_y3_ct2ved__4(inter->w_phys_ct, scratch2, par);
	
	for (i=0; i<=local_size_ved-1; ++i){
		scratch1[i] = ( 9.0*scratch1[i] - 1.0*scratch2[i] )/8.0;
	}
	
	
    // ved is added at j = 0, ny
	// but from i=start to end only.
	// needs i=start-2,start-1,end+1
    get_dwdy_ved_wall_ver2(inter->w_phys_ct, scratch1, par);
	
	interpolate_x1_ved2cn__4(scratch1, var->dwdy, par);
	interpolate_x3_ved2cn__4(scratch1, scratch3, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->dwdy[i] = ( 9.0*var->dwdy[i] - 1.0*scratch3[i] )/8.0;
	}
	
	
	fftw_free(scratch1);
	fftw_free(scratch2);
	fftw_free(scratch3);
	
	
	
	
	// get dudy *******************************************
	// need some supplimental work at j = 0, ny-1
	// since dwdy at j = 0, ny-1 are set to be zero
	// using the following methods.
	scratch1 = FFTW_MALLOC(local_size_cn);
	derivative_y1_ued2cn__4(inter->u_phys_ued, var->dudy, par);
	derivative_y3_ued2cn__4(inter->u_phys_ued, scratch1, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->dudy[i] = ( 9.0*var->dudy[i] - 1.0*scratch1[i] )/8.0;
	}
	
	get_dudy_cn_wall(inter->u_phys_ued, u_t, u_b, var->dudy, par);
	
	// get dvdx *******************************************
	derivative_x1_ved2cn__4(inter->v_phys_ved, var->dvdx, par);
	derivative_x3_ved2cn__4(inter->v_phys_ved, scratch1, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		var->dvdx[i] = ( 9.0*var->dvdx[i] - 1.0*scratch1[i] )/8.0;
	}
	
	fftw_free(scratch1);
	
	fftw_free(u_t);
	fftw_free(u_b);
	fftw_free(w_t);
	fftw_free(w_b);
	
	
	
	// At the wall, the value should be zero
	
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			var->dvdx[par->index_cn[i][0][k]] = 0.0;
			var->dvdx[par->index_cn[i][ny][k]] = 0.0;
			
			var->dvdz[par->index_cn[i][0][k]] = 0.0;
			var->dvdz[par->index_cn[i][ny][k]] = 0.0;
			
			var->dwdx[par->index_cn[i][0][k]] = 0.0;
			var->dwdx[par->index_cn[i][ny][k]] = 0.0;
			
			var->dwdz[par->index_cn[i][0][k]] = 0.0;
			var->dwdz[par->index_cn[i][ny][k]] = 0.0;
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/* get velocity gradient tensor defined at corner					   */
/* these values will be used to calculate subgrid stress tensor		   */
/* Originally,														   */
/* dudx, dvdy and dwdz are defined at center						   */
/* will be interpolated to v-edge, then to corner					   */
/* dudz, dwdx are defined at u-edge, will be interpolated to corner    */
/* dvdz, dwdy are defined at v-edge, will be interpolated to corner    */
/* dudy, dvdx are defined at corner (no interpolation needed)		   */
/* ******************************************************************* */
void get_velgrad_tensor_cn__2(
                              variables *var,
                              interpolated *inter,
                              parameters *par,
                              fftwplans *planptr)
{
	
	int i, k;
	
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	double *scratch1;
	
	
	scratch1 = FFTW_MALLOC(local_size_ct);
	
	
	// get dudx *******************************************
    // ued[start-1:end][-1:ny] --> ct[start-1:end-1][-1:ny]
    // ct[end] is not defined
	derivative_x1_ued2ct__2(inter->u_phys_ued, scratch1, par);
	
    // ct[start-1:end-1][-1:ny] --> cn[start:end-1][0:ny]
	interpolate_ct2cn__2(scratch1, var->dudx, par);
	
	
	// get dvdy *******************************************
	// values at j = -1, ny are not defined in scratch1
	derivative_y1_ved2ct__2(inter->v_phys_ved, scratch1, par);
	
	interpolate_ct2cn__2(scratch1, var->dvdy, par);
	
	// in addition, fix dvdy at the wall
	{
		
		double *dvdy_ved;
		double *scratch_cn1;
		dvdy_ved = FFTW_MALLOC(local_size_ved);
		scratch_cn1 = FFTW_MALLOC(local_size_cn);
		
		for (i=0; i<=local_size_ved-1; ++i){
			dvdy_ved[i] = 0.0;
		}
		
		for (i=start-1; i<=end; ++i){
			for (k=0; k<=nz-1; ++k){
				dvdy_ved[par->index_ved[i][0][k]]
				= (-   11.0*inter->v_phys_ved[par->index_ved[i][0][k]]
				   +   18.0*inter->v_phys_ved[par->index_ved[i][1][k]]
				   -    9.0*inter->v_phys_ved[par->index_ved[i][2][k]]
				   +    2.0*inter->v_phys_ved[par->index_ved[i][3][k]])/(6.0*par->dy);
				
				dvdy_ved[par->index_ved[i][ny][k]]
				= (+   11.0*inter->v_phys_ved[par->index_ved[i][ny][k]]
				   -   18.0*inter->v_phys_ved[par->index_ved[i][ny-1][k]]
				   +    9.0*inter->v_phys_ved[par->index_ved[i][ny-2][k]]
				   -    2.0*inter->v_phys_ved[par->index_ved[i][ny-3][k]])/(6.0*par->dy);
				
			}
		}
		
		interpolate1_ved_cn__2(dvdy_ved, scratch_cn1, par);
		
		for (i=start; i<=end; ++i){
			for (k=0; k<=nz-1; ++k){
				var->dvdy[par->index_cn[i][0][k]]
				= scratch_cn1[par->index_cn[i][0][k]];
				
				var->dvdy[par->index_cn[i][ny][k]]
				= scratch_cn1[par->index_cn[i][ny][k]];
				
			}
		}
		
		
		fftw_free(dvdy_ved);
		fftw_free(scratch_cn1);
		
	}
	
	// get dwdz *******************************************
	derivative_z_ct2ct__2(inter->w_ct, scratch1, par);
	fftw_execute_r2r(planptr->p1d_invz_ct, scratch1, scratch1);
	
	interpolate_ct2cn__2(scratch1, var->dwdz, par);
	
    fftw_free(scratch1);
	
	
	// get dudz ******************************************
	scratch1 = FFTW_MALLOC(local_size_ued);
	derivative_z_ued2ued__2(inter->u_ued, scratch1, par);
	
	fftw_execute_r2r(planptr->p1d_invz_ued, scratch1, scratch1);
	
	interpolate1_ued_cn__2(scratch1, var->dudz, par);
	
	
	
	// get dwdx ******************************************
	derivative_x1_ct2ued__2(inter->w_phys_ct, scratch1, par);
	
	interpolate1_ued_cn__2(scratch1, var->dwdx, par);
	
	
	
	fftw_free(scratch1);
	
	
	
	
	
	// get dvdz ******************************************
	scratch1 = FFTW_MALLOC(local_size_ved);
	derivative_z_ved2ved__2(inter->v_ved, scratch1, par);
	
	fftw_execute_r2r(planptr->p1d_invz_ved, scratch1, scratch1);
	
	interpolate1_ved_cn__2(scratch1, var->dvdz, par);
	
	
	// get dwdy ******************************************
	derivative_y1_ct2ved__2(inter->w_phys_ct, scratch1, par);
	interpolate1_ved_cn__2(scratch1, var->dwdy, par);
	
	
	fftw_free(scratch1);
	
	// get dudy *******************************************
	derivative_y1_ued2cn__2(inter->u_phys_ued, var->dudy, par);
	
	
	// get dvdx *******************************************
	derivative_x1_ved2cn__2(inter->v_phys_ved, var->dvdx, par);
	
    
	
	return;
}






void get_dudy_cn_wall(double *u, double *u_t, double *u_b, double *dudy, parameters *par)
{
	
	int i, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start+ par->local_nx;
	
	double u_top, u_bot;
	
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			u_top = 3.0/8.0*u[par->index_ued[i][ny][k]]
            + 3.0/4.0*u[par->index_ued[i][ny-1][k]]
            - 1.0/8.0*u[par->index_ued[i][ny-2][k]];
			
			dudy[par->index_cn[i][ny][k]] \
            = (+ 184*u_top
               - 225*u[par->index_ued[i][ny-1][k]]
               +  50*u[par->index_ued[i][ny-2][k]]
               -   9*u[par->index_ued[i][ny-3][k]] )/(60.0*par->dy);
            
			u_bot = 3.0/8.0*u[par->index_ued[i][-1][k]]
            + 3.0/4.0*u[par->index_ued[i][0][k]]
            - 1.0/8.0*u[par->index_ued[i][1][k]];
			
			dudy[par->index_cn[i][0][k]] \
            = (- 184*u_bot
               + 225*u[par->index_ued[i][0][k]]
               -  50*u[par->index_ued[i][1][k]]
               +   9*u[par->index_ued[i][2][k]])/(60.0*par->dy);
    	}
	}
	
	return;
}


void get_dudy_wp_wall(double *u, double *u_t, double *u_b, double *dudy, parameters *par)
{
	
	int i, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double u_top, u_bot;
	double *dudy_ued;
	double *scratch;
	
	scratch = FFTW_MALLOC(par->local_size_wp);
	dudy_ued = FFTW_MALLOC(par->local_size_ued);
	
	for (i=0; i<=par->local_size_ued-1; ++i){
		dudy_ued[i] = 0.0;
	}
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			u_top = 3.0/8.0*u[par->index_ued[i][ny][k]]
            + 3.0/4.0*u[par->index_ued[i][ny-1][k]]
            - 1.0/8.0*u[par->index_ued[i][ny-2][k]];
            
			dudy_ued[par->index_ued[i][ny-1][k]] \
            = (+ 32*u_top
               - 15*u[par->index_ued[i][ny-1][k]]
               - 20*u[par->index_ued[i][ny-2][k]]
               +  3*u[par->index_ued[i][ny-3][k]] )/(30.0*par->dy);
			
			
			u_bot = 3.0/8.0*u[par->index_ued[i][-1][k]]
            + 3.0/4.0*u[par->index_ued[i][0][k]]
            - 1.0/8.0*u[par->index_ued[i][1][k]];
			
			dudy_ued[par->index_ued[i][0][k]] \
            = (- 32*u_bot
               + 15*u[par->index_ued[i][0][k]]
               + 20*u[par->index_ued[i][1][k]]
               -  3*u[par->index_ued[i][2][k]])/(30.0*par->dy);
		}
	}
	
	interpolate_x1_ued2wp__4(dudy_ued, scratch, par);
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			dudy[par->index_wp[i][0][k]]
            = 9.0/8.0*scratch[par->index_wp[i][0][k]];
			dudy[par->index_wp[i][ny-1][k]]
            = 9.0/8.0*scratch[par->index_wp[i][ny-1][k]];
		}
	}
	
	interpolate_x3_ued2wp__4(dudy_ued, scratch, par);
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			dudy[par->index_wp[i][0][k]]
			+= - 1.0/8.0*scratch[par->index_wp[i][0][k]];
			dudy[par->index_wp[i][ny-1][k]]
			+= - 1.0/8.0*scratch[par->index_wp[i][ny-1][k]];
		}
	}
	
	
	
	
	return;
}


void get_dwdy_ved_wall_ver2(double *w, double *dwdy, parameters *par)
{
	
	int i, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
    const double inv60dy = 1.0/(60*par->dy);
    for (i=start-2; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
            
            dwdy[par->index_ved[i][0][k]] \
            = (- 69*w[par->index_ct[i][-1][k]]
               + 87*w[par->index_ct[i][0][k]]
               - 27*w[par->index_ct[i][1][k]]
               +  9*w[par->index_ct[i][2][k]])*inv60dy;
            
            
            dwdy[par->index_ved[i][ny][k]] \
            = (+ 69*w[par->index_ct[i][ny][k]]
               - 87*w[par->index_ct[i][ny-1][k]]
               + 27*w[par->index_ct[i][ny-2][k]]
               -  9*w[par->index_ct[i][ny-3][k]])*inv60dy;
            
			
		}
	}
	return;
}


void get_dwdy_wp_wall(double *w, double *w_t, double *w_b, double *dwdy, parameters *par)
{
	
	int i, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			dwdy[par->index_wp[i][ny-1][k]] \
            = (+ 32.0*w_t[par->index_tb[i][k]]
               - 15.0*w[par->index_ct[i][ny-1][k]]
               - 20.0*w[par->index_ct[i][ny-2][k]]
               +  3.0*w[par->index_ct[i][ny-3][k]] )/(30.0*par->dy);
			
			dwdy[par->index_wp[i][0][k]] \
            = (- 32.0*w_b[par->index_tb[i][k]]
               + 15.0*w[par->index_ct[i][0][k]]
               + 20.0*w[par->index_ct[i][1][k]]
               -  3.0*w[par->index_ct[i][2][k]])/(30.0*par->dy);
		}
	}
	
	
	return;
}








