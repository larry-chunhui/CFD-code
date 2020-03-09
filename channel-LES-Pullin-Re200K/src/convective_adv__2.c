/*
 *  convective_adv__2.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Calculates the advective form of the convective terms with
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
#include "differentiate.h"
#include "convective.h"
#include "convective_div__2.h"
#include "convective_adv__2.h"

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
void get_convective_adv__2(
                           interpolated *inter,
                           double *Nu_copy,
                           double *Nv_copy,
                           double *Nw_copy,
                           fftwplans *planptr,
                           parameters *par)
{
	int i, j, k;
	int index;
	int index1, index2, index3, index4;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	double *vdudy, *udvdx; // defined at index_cn
	double *udwdx, *wdudz; // defined at index_ued
	double *wdvdz, *vdwdy; // defined at index_ved
	double *vdvdy, *ududx, *wdwdz; // defined at index_wp
	
	// for LES wall ODE
	double *wall_ududx = NULL, *wall_wdudz = NULL;
	if (par->les == 1 && par->bc == 1){
		wall_ududx = FFTW_MALLOC(par->local_size_u);
		wall_wdudz = FFTW_MALLOC(par->local_size_u);
	}
    
	
	vdudy = FFTW_MALLOC(local_size_cn); // for Nv
	udvdx = FFTW_MALLOC(local_size_cn); // for Nu
	
	udwdx = FFTW_MALLOC(local_size_ued); // for Nw
	wdudz = FFTW_MALLOC(local_size_ued); // for Nu
	
	wdvdz = FFTW_MALLOC(local_size_ved); // for Nv
	vdwdy = FFTW_MALLOC(local_size_ved); // for Nw
	
	vdvdy = FFTW_MALLOC(local_size_ct); // for Nv
	ududx = FFTW_MALLOC(local_size_ct); // for Nx
	wdwdz = FFTW_MALLOC(local_size_ct); // for Nz
    
	for (i=0; i<=local_size_cn-1; ++i){
		vdudy[i] = 0.0;
		udvdx[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ued-1; ++i){
		udwdx[i] = 0.0;
		wdudz[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ved-1; ++i){
		wdvdz[i] = 0.0;
		vdwdy[i] = 0.0;
	}
	
	for (i=0; i<=local_size_ct-1; ++i){
		ududx[i] = 0.0;
		vdvdy[i] = 0.0;
		wdwdz[i] = 0.0;
	}
	
	form_product_adi(inter, \
                     ududx, udvdx, udwdx, \
                     vdudy, vdvdy, vdwdy, \
                     wdudz, wdvdz, wdwdz, \
                     planptr, par);
    
	
	// Nu
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Nu_copy[par->index_u[i][j][k]] \
                = (ududx[par->index_ct[i-1][j][k]] + ududx[par->index_ct[i][j][k]])/2.0 \
                + (vdudy[par->index_cn[i][j][k]] + vdudy[par->index_cn[i][j+1][k]])/2.0 \
                +  wdudz[par->index_ued[i][j][k]];
				
				// for LES wall ODE
				if (par->les == 1 && par->bc == 1){
					wall_ududx[par->index_u[i][j][k]]
					= (ududx[par->index_ct[i-1][j][k]] + ududx[par->index_ct[i][j][k]])/2.0;
					wall_wdudz[par->index_u[i][j][k]]
					=  wdudz[par->index_ued[i][j][k]];
				}
				
			}
		}
	}
	
	// Nv
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Nv_copy[par->index_v[i][j][k]] \
                = (udvdx[par->index_cn[i][j][k]] + udvdx[par->index_cn[i+1][j][k]])/2.0 \
                + (vdvdy[par->index_ct[i][j-1][k]] + vdvdy[par->index_ct[i][j][k]])/2.0 \
                +  wdvdz[par->index_ved[i][j][k]];
			}
		}
	}
	
	// Nw
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Nw_copy[par->index_wp[i][j][k]] \
                = (udwdx[par->index_ued[i][j][k]] + udwdx[par->index_ued[i+1][j][k]])/2.0 \
                + (vdwdy[par->index_ved[i][j][k]] + vdwdy[par->index_ved[i][j+1][k]])/2.0 \
                +  wdwdz[par->index_ct[i][j][k]];
			}
		}
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
					inter->bc_ududx_b[index] = wall_ududx[index1];
					inter->bc_ududx_t[index] = wall_ududx[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_ududx_b[index]
					= (wall_ududx[index1]+wall_ududx[index3])/2.0;
					inter->bc_ududx_t[index]
					= (wall_ududx[index2]+wall_ududx[index4])/2.0;
				}
			}
		}
	}
	// *****************************************************************
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
					inter->bc_wdudz_b[index] = wall_wdudz[index1];
					inter->bc_wdudz_t[index] = wall_wdudz[index2];
				}
				
				if (par->ode_choice == 2) {
					inter->bc_wdudz_b[index]
					= (wall_wdudz[index1]+wall_wdudz[index3])/2.0;
					inter->bc_wdudz_t[index]
					= (wall_wdudz[index2]+wall_wdudz[index4])/2.0;
				}
			}
		}
		
		
		// transfer into fourier space
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_wdudz_b, inter->bc_wdudz_b);
		fftw_execute_r2r(planptr->p1d_z_bc_tb, inter->bc_wdudz_t, inter->bc_wdudz_t);
		for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
			inter->bc_wdudz_t[i] = inter->bc_wdudz_t[i]/(double)(nz);
			inter->bc_wdudz_b[i] = inter->bc_wdudz_b[i]/(double)(nz);
		}
	}
	// *****************************************************************
	
	if (par->les == 1 && par->bc == 1){
		fftw_free(wall_ududx);
		fftw_free(wall_wdudz);
	}
	
	
	// free allocated memory
	fftw_free(ududx);
	fftw_free(udvdx);
	fftw_free(udwdx);
	fftw_free(vdudy);
	fftw_free(vdvdy);
	fftw_free(vdwdy);
	fftw_free(wdudz);
	fftw_free(wdvdz);
	fftw_free(wdwdz);
	
	return;
}

/* ******************************************************************* */

/* ******************************************************************* */
void form_product_adi(
                      interpolated *inter,
                      double *ududx,
                      double *udvdx,
                      double *udwdx,
                      double *vdudy,
                      double *vdvdy,
                      double *vdwdy,
                      double *wdudz,
                      double *wdvdz,
                      double *wdwdz,
                      fftwplans *planptr,
                      parameters *par)
{
	
	form_product_x(inter, ududx, udvdx, udwdx, planptr, par);
	form_product_y(inter, vdudy, vdvdy, vdwdy, planptr, par);
	form_product_z(inter, wdudz, wdvdz, wdwdz, planptr, par);
	
	return;
}

/* ******************************************************************* */

/* ******************************************************************* */
void form_product_x(
                    interpolated *inter,
                    double *ududx,
                    double *udvdx,
                    double *udwdx,
                    fftwplans *planptr,
                    parameters *par)
{
    
	int i;
	const int nz = par->nz;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
    
	double *dudx;
	double *dvdx;
	double *dwdx;
    
	dudx = FFTW_MALLOC(local_size_ct);
	dvdx = FFTW_MALLOC(local_size_cn);
	dwdx = FFTW_MALLOC(local_size_ued);
	
	
	get_dudx2(inter->u_phys_ued, dudx, par); // dudx defined at ct
	get_dvdx2(inter->v_phys_ved, dvdx, par); // dvdx defined at cn
	get_dwdx2(inter->w_phys_ct,  dwdx, par); // dwdx defined at ued
	
	calc_multiply_ct( dudx,	inter->u_phys_ct1, ududx, par);
	calc_multiply_cn( dvdx,	inter->u_phys_cn1, udvdx, par);
	calc_multiply_ued(dwdx, inter->u_phys_ued, udwdx, par);
	
	fftw_execute_r2r(planptr->p1d_z_ct,  ududx, ududx);
	fftw_execute_r2r(planptr->p1d_z_cn,  udvdx, udvdx);
	fftw_execute_r2r(planptr->p1d_z_ued, udwdx, udwdx);
    
	for (i=0; i<=local_size_ct-1; ++i){
		ududx[i] = ududx[i]/(double)nz;
	}
    
	for (i=0; i<=local_size_cn-1; ++i){
		udvdx[i] = udvdx[i]/(double)nz;
	}
	
	for (i=0; i<=local_size_ued-1; ++i){
		udwdx[i] = udwdx[i]/(double)nz;
	}
	
	
	fftw_free(dudx);
	fftw_free(dvdx);
	fftw_free(dwdx);
	
	
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void form_product_y(
                    interpolated *inter,
                    double *vdudy,
                    double *vdvdy,
                    double *vdwdy,
                    fftwplans *planptr,
                    parameters *par)
{
    
    
	int i;
	const int nz = par->nz;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
    
	double *dudy;
	double *dvdy;
	double *dwdy;
    
	dudy = FFTW_MALLOC(local_size_cn);
	dvdy = FFTW_MALLOC(local_size_ct);
	dwdy = FFTW_MALLOC(local_size_ved);
	
	get_dudy2(inter->u_phys_ued, dudy, par); // dudx defined at cn
	get_dvdy2(inter->v_phys_ved, dvdy, par); // dvdx defined at ct
	get_dwdy2(inter->w_phys_ct,  dwdy, par); // dwdx defined at ved
	
	calc_multiply_cn( dudy,	inter->v_phys_cn1,  vdudy, par);
	calc_multiply_ct( dvdy,	inter->v_phys_ct1,  vdvdy, par);
	calc_multiply_ved(dwdy, inter->v_phys_ved, vdwdy, par);
	
	fftw_execute_r2r(planptr->p1d_z_cn,  vdudy, vdudy);
	fftw_execute_r2r(planptr->p1d_z_ct,  vdvdy, vdvdy);
	fftw_execute_r2r(planptr->p1d_z_ved, vdwdy, vdwdy);
    
	for (i=0; i<=local_size_cn-1; ++i){
		vdudy[i] = vdudy[i]/(double)nz;
	}
    
	for (i=0; i<=local_size_ct-1; ++i){
		vdvdy[i] = vdvdy[i]/(double)nz;
	}
	
	for (i=0; i<=local_size_ved-1; ++i){
		vdwdy[i] = vdwdy[i]/(double)nz;
	}
	
	
	fftw_free(dudy);
	fftw_free(dvdy);
	fftw_free(dwdy);
	
	return;
}


/* ******************************************************************* */

/* ******************************************************************* */
void form_product_z(
                    interpolated *inter,
                    double *wdudz,
                    double *wdvdz,
                    double *wdwdz,
                    fftwplans *planptr,
                    parameters *par)
{
	
	int i;
	const int nz = par->nz;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	
	double *dudz;
	double *dvdz;
	double *dwdz;
    
	dudz = FFTW_MALLOC(local_size_ued);
	dvdz = FFTW_MALLOC(local_size_ved);
	dwdz = FFTW_MALLOC(local_size_ct);
	
	get_dudz_fft(inter->u_ued, dudz, par);
	get_dvdz_fft(inter->v_ved, dvdz, par);
	get_dwdz_fft(inter->w_ct,  dwdz, par);
    
	fftw_execute_r2r(planptr->p1d_invz_ued, dudz, dudz);
	fftw_execute_r2r(planptr->p1d_invz_ved, dvdz, dvdz);
	fftw_execute_r2r(planptr->p1d_invz_ct,  dwdz, dwdz);
	
	calc_multiply_ued(dudz, inter->w_phys_ued1, wdudz, par);
	calc_multiply_ved(dvdz, inter->w_phys_ved1, wdvdz, par);
	calc_multiply_ct( dwdz, inter->w_phys_ct,  wdwdz, par);
	
	fftw_execute_r2r(planptr->p1d_z_ued, wdudz, wdudz);
	fftw_execute_r2r(planptr->p1d_z_ved, wdvdz, wdvdz);
	fftw_execute_r2r(planptr->p1d_z_ct,  wdwdz, wdwdz);
    
	for (i=0; i<=local_size_ued-1; ++i){
		wdudz[i] = wdudz[i]/(double)nz;
	}
    
	for (i=0; i<=local_size_ved-1; ++i){
		wdvdz[i] = wdvdz[i]/(double)nz;
	}
	
	for (i=0; i<=local_size_ct-1; ++i){
		wdwdz[i] = wdwdz[i]/(double)nz;
	}
	
	
	fftw_free(dudz);
	fftw_free(dvdz);
	fftw_free(dwdz);
	
	return;
}




/* ******************************************************************* */
/* take the derivative of u(index_ued) w.r.t. x (streamwise direction) */
/* derivatives are defined at center, index_ct						   */
/* dudx is not defined at i = end									   */
/* ******************************************************************* */
void get_dudx2(
               double *u,
               double *dudx,
               parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double dx = par->dx;
	
	//interior
	for (i=start-1; i<=end-1; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dudx[par->index_ct[i][j][k]] \
                = (u[par->index_ued[i+1][j][k]] - u[par->index_ued[i][j][k]])/dx;
			}
		}
	}
	
	i = end;
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			dudx[par->index_ct[i][j][k]] = 0.0;
		}
	}
    
    
	return;
}


/* ******************************************************************* */
/* take the derivative of v(index_ved) w.r.t. x (streamwise direction) */
/* derivatives are defined at corner (index_cn)						   */
/* ******************************************************************* */
void get_dvdx2(
               double *v,
               double *dvdx,
               parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double dx = par->dx;
	
	for (i=start; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdx[par->index_cn[i][j][k]] \
                = (v[par->index_ved[i][j][k]] - v[par->index_ved[i-1][j][k]])/dx;
			}
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/* take the derivative of w(index_ct) w.r.t. x (streamwise direction)  */
/* derivatives are defined at index_ued								   */
/* dwdx is not defined at i = start-1								   */
/* ******************************************************************* */
void get_dwdx2(
               double *w,
               double *dwdx,
               parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double dx = par->dx;
	
	for (i=start; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdx[par->index_ued[i][j][k]] \
                = (w[par->index_ct[i][j][k]] - w[par->index_ct[i-1][j][k]])/dx;
			}
		}
	}
	
	i = start-1;
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			dwdx[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	
    
	return;
}



/* ******************************************************************* */
/* take the derivative of u(index_ued) w.r.t. y (wall-normal direction)*/
/* derivatives are defined at corner (index_cn)						   */
/* ******************************************************************* */
void get_dudy2(
               double *u,
               double *dudy,
               parameters *par)
{
    
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double dy = par->dy;
	
	for (i=start; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dudy[par->index_cn[i][j][k]] \
                = (u[par->index_ued[i][j][k]] - u[par->index_ued[i][j-1][k]])/dy;
			}
		}
	}
    
	
	return;
}


/* ******************************************************************* */
/* take the derivative of v(index_ved) w.r.t. y (wall-normal direction)*/
/* derivatives are defined at center, index_ct  					   */
/* dvdy is not defined at j = -1, ny								   */
/* ******************************************************************* */
void get_dvdy2(
               double *v,
               double *dvdy,
               parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double dy = par->dy;
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdy[par->index_ct[i][j][k]] \
                = (v[par->index_ved[i][j+1][k]] - v[par->index_ved[i][j][k]])/dy;
			}
		}
	}
	
	j = -1;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			dvdy[par->index_ct[i][j][k]] = 0.0;
		}
	}
	
	j = ny;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			dvdy[par->index_ct[i][j][k]] = 0.0;
		}
	}
	
	return;
}


/* ******************************************************************* */
/* take the derivative of w(index_ct) w.r.t. y (wall-normal direction) */
/* derivatives are defined at v edge, index_ved						   */
/* ******************************************************************* */
void get_dwdy2(
               double *w,
               double *dwdy,
               parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double dy = par->dy;
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdy[par->index_ved[i][j][k]] \
                = (w[par->index_ct[i][j][k]] - w[par->index_ct[i][j-1][k]])/dy;
			}
		}
	}
	
	return;
}



/* ******************************************************************* */
/* take the derivative of u(index_ued) w.r.t. w (fft direction)        */
/* derivatives are also defined at u edge, index_ued				   */
/* ******************************************************************* */
void get_dudz_fft(
                  double *u,
                  double *dudz,
                  parameters *par)
{
    
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	for (i=start-1; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=1; k<=(nz-1)/2; ++k){
				dudz[par->index_ued[i][j][nz-k]] = u[par->index_ued[i][j][k]]*par->kz[k];
				dudz[par->index_ued[i][j][k]] =  - u[par->index_ued[i][j][nz-k]]*par->kz[nz-k];
			}
			k = 0;
			dudz[par->index_ued[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				dudz[par->index_ued[i][j][k]] = 0.0;
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/* take the derivative of v (index_ved) w.r.t. w (fft direction)       */
/* derivatives are also defined at v edge, index_ved				   */
/* ******************************************************************* */
void get_dvdz_fft(
                  double *v,
                  double *dvdz,
                  parameters *par)
{
    
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			for (k=1; k<=(nz-1)/2; ++k){
				dvdz[par->index_ved[i][j][nz-k]] = v[par->index_ved[i][j][k]]*par->kz[k];
				dvdz[par->index_ved[i][j][k]] =  - v[par->index_ved[i][j][nz-k]]*par->kz[nz-k];
			}
			k = 0;
			dvdz[par->index_ved[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				dvdz[par->index_ved[i][j][k]] = 0.0;
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/* take the derivative of w (defined at ct) w.r.t. w (fft direction)   */
/* derivatives are also defined at center, index_ct					   */
/* ******************************************************************* */
void get_dwdz_fft(
                  double *w,
                  double *dwdz,
                  parameters *par)
{
    
	int i, j, k;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	for (i=start-1; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=1; k<=(nz-1)/2; ++k){
				dwdz[par->index_ct[i][j][nz-k]] = w[par->index_ct[i][j][k]]*par->kz[k];
				dwdz[par->index_ct[i][j][k]] =  - w[par->index_ct[i][j][nz-k]]*par->kz[nz-k];
			}
			k = 0;
			dwdz[par->index_ct[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				dwdz[par->index_ct[i][j][k]] = 0.0;
			}
		}
	}
	
	return;
}

