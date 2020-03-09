/*
 *  interpolation__4.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Fourth order interpolation.
 */

#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "interpolation.h"
#include "interpolation__4.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


// note: can be integrated with 2nd order subroutines
void interpolation__4(
                      int flag,
                      variables *var,
                      neighbors *shared,
                      interpolated *inter,
                      parameters *par,
                      fftwplans *planptr)
{
	
	int i, j, k;
	int index1, index2;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	double *copy_ct;
	double *copy_ued;
	double *copy_ved;
	
	bigger_array_u_ued__4(flag, var, shared, inter->u_ued, par); // checked
	bigger_array_v_ved__4(flag, var, shared, inter->v_ved, par); // checked
	bigger_array_wp_ct__4(flag, var, shared, inter->w_ct, par); // checked
	
	// obtained from continuity relation at cells j = -1
	bigger_array_v_ved__4_addition(inter->u_ued, inter->v_ved, inter->w_ct, par);
	
	
	//************************* index_ct ***************************
	copy_ct = FFTW_MALLOC(par->local_size_ct);
	
	for (i=0; i<=par->local_size_ct-1; ++i){
		copy_ct[i] = inter->w_ct[i];
	}
	fftw_execute_r2r(planptr->p1d_invz_ct, copy_ct, inter->w_phys_ct);
	
	fftw_free(copy_ct);
	
	//************************* index_ued ***************************
	copy_ued = FFTW_MALLOC(par->local_size_ued);
	
	for (i=0; i<=par->local_size_ued-1; ++i){
		copy_ued[i] = inter->u_ued[i];
	}
	fftw_execute_r2r(planptr->p1d_invz_ued, copy_ued, inter->u_phys_ued);
	
	fftw_free(copy_ued);
	
	//************************* index_ved ***************************
	copy_ved = FFTW_MALLOC(par->local_size_ved);
	
	for (i=0; i<=par->local_size_ved-1; ++i){
		copy_ved[i] = inter->v_ved[i];
	}
	fftw_execute_r2r(planptr->p1d_invz_ved, copy_ved, inter->v_phys_ved);
	
	fftw_free(copy_ved);
	
	
	//************************* interpolation ***************************
	interpolate_y1_ued2cn__4(inter->u_phys_ued, inter->u_phys_cn1, par);
	interpolate_x1_ued2ct__4(inter->u_phys_ued, inter->u_phys_ct1, par);
	
	interpolate_x1_ved2cn__4(inter->v_phys_ved, inter->v_phys_cn1, par);
	interpolate_y1_ved2ct__4(inter->v_phys_ved, inter->v_phys_ct1, par);
	
	interpolate_x1_ct2ued__4(inter->w_phys_ct, inter->w_phys_ued1, par);
	interpolate_y1_ct2ved__4(inter->w_phys_ct, inter->w_phys_ved1, par);
	
	//************************* interpolation ***************************
	interpolate_y3_ued2cn__4(inter->u_phys_ued, inter->u_phys_cn3, par);
	interpolate_x3_ued2ct__4(inter->u_phys_ued, inter->u_phys_ct3, par);
	
	interpolate_x3_ved2cn__4(inter->v_phys_ved, inter->v_phys_cn3, par);
	interpolate_y3_ved2ct__4(inter->v_phys_ved, inter->v_phys_ct3, par);
	
	interpolate_x3_ct2ued__4(inter->w_phys_ct, inter->w_phys_ued3, par);
	interpolate_y3_ct2ved__4(inter->w_phys_ct, inter->w_phys_ved3, par);
	
	
	
	double *w_phys_cn;
	
	if (par->les == 1){
        
		w_phys_cn = FFTW_MALLOC(par->local_size_cn);
        
		interpolate_ct2cn__4(inter->w_phys_ct, w_phys_cn, par);
        
		for (i=start-1; i<=end; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					index1 = par->index_les[i][j][k];
					index2 = par->index_cn[i][j][k];
                    
					inter->u_phys_les[index1]
					= (9.0*inter->u_phys_cn1[index2] - 1.0*inter->u_phys_cn3[index2])/8.0;
                    
					inter->v_phys_les[index1]
					= (9.0*inter->v_phys_cn1[index2] - 1.0*inter->v_phys_cn3[index2])/8.0;
                    
					inter->w_phys_les[index1] = w_phys_cn[index2];
					
				}
			}
		}
		
        
		fftw_free(w_phys_cn);
	}
	
    
	return;
}



/* ******************************************************************* */
/*
 interpolate velocity w[] (j & j-1 values for j)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 
 w_ved is not defined
 at j = -2, ny+2
 */
/* ******************************************************************* */
void interpolate_y1_ct2ved__4(
                              double *w_ct,
                              double *w_ved,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-3; i<=end+2; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_ved[par->index_ved[i][j][k]] \
                = (w_ct[par->index_ct[i][j-1][k]] + w_ct[par->index_ct[i][j][k]])/2.0;
			}
		}
	}
	
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			w_ved[par->index_ved[i][-2][k]] = 0.0;
			w_ved[par->index_ved[i][ny+2][k]] = 0.0;
		}
	}
	
	
	
	return;
}




/* ******************************************************************* */
/*
 interpolate velocity w[] (j+1 & j-2 values for j)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 
 w_ved is not defined
 at j = -2, -1, ny+1, ny+2
 */
/* ******************************************************************* */
void interpolate_y3_ct2ved__4(
                              double *w_ct,
                              double *w_ved,
                              parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-3; i<=end+2; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				w_ved[par->index_ved[i][j][k]] \
                = (w_ct[par->index_ct[i][j-2][k]] + w_ct[par->index_ct[i][j+1][k]])/2.0;
			}
		}
	}
	
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			//j = -2;
			w_ved[par->index_ved[i][-2][k]] = 0.0;
			//j = -1;
			w_ved[par->index_ved[i][-1][k]] = 0.0;
		}
	}
	
    
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			//j = ny+1;
			w_ved[par->index_ved[i][ny+1][k]] = 0.0;
			//j = ny+2;
			w_ved[par->index_ved[i][ny+2][k]] = 0.0;
		}
	}
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity w[] (i-1 & i values for i)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 w_ued is not defined
 at i = start-3
 */
/* ******************************************************************* */
void interpolate_x1_ct2ued__4(
                              double *w_ct,
                              double *w_ued,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-2; i<=end+2; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_ued[par->index_ued[i][j][k]] \
                = (w_ct[par->index_ct[i-1][j][k]] + w_ct[par->index_ct[i][j][k]])/2.0;
			}
		}
	}
	
	
	i = start-3;
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			w_ued[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity w[] (i-2 & i+1 values for i)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 w_ued is not defined
 at i = start-3, start-2, end+2
 */
/* ******************************************************************* */
void interpolate_x3_ct2ued__4(
                              double *w_ct,
                              double *w_ued,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_ued[par->index_ued[i][j][k]] \
                = (w_ct[par->index_ct[i-2][j][k]] + w_ct[par->index_ct[i+1][j][k]])/2.0;
			}
		}
	}
	
	
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			//i = start-3;
			w_ued[par->index_ued[start-3][j][k]] = 0.0;
			//i = start-2;
			w_ued[par->index_ued[start-2][j][k]] = 0.0;
			//i = end+2;
			w_ued[par->index_ued[end+2][j][k]] = 0.0;
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity v[] (i-1 & i values for i)
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_x1_ved2cn__4(
                              double *v_ved,
                              double *v_cn,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				v_cn[par->index_cn[i][j][k]] \
                = (v_ved[par->index_ved[i-1][j][k]] + v_ved[par->index_ved[i][j][k]])/2.0;
			}
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/*
 interpolate velocity v[] (i-2 & i+1 values for i)
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_x3_ved2cn__4(
                              double *v_ved,
                              double *v_cn,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				v_cn[par->index_cn[i][j][k]] \
                = (v_ved[par->index_ved[i-2][j][k]] + v_ved[par->index_ved[i+1][j][k]])/2.0;
			}
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity v[] (index_v[][][]) (j & j+1 values for j)
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_y1_ved2ct__4(
                              double *v_ved,
                              double *v_ct,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start-3; i<=end+2; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				v_ct[par->index_ct[i][j][k]] \
                = (v_ved[par->index_ved[i][j][k]] + v_ved[par->index_ved[i][j+1][k]])/2.0;
			}
		}
	}
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity v[] (index_v[][][]) (j-1 & j+2 values for j)
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 v_ct is not defined at j = -2, j = ny+1
 */
/* ******************************************************************* */
void interpolate_y3_ved2ct__4(
                              double *v_ved,
                              double *v_ct,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start-3; i<=end+2; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				v_ct[par->index_ct[i][j][k]] \
                = (v_ved[par->index_ved[i][j-1][k]] + v_ved[par->index_ved[i][j+2][k]])/2.0;
			}
		}
	}
	
	//j = -2;
	//j = ny+1;
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			v_ct[par->index_ct[i][-2][k]] = 0.0;
			v_ct[par->index_ct[i][ny+1][k]] = 0.0;
		}
	}
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity u[] (j-1 & j values for j)
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_y1_ued2cn__4(
                              double *u_ued,
                              double *u_cn,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_cn[par->index_cn[i][j][k]] \
                = (u_ued[par->index_ued[i][j-1][k]] + u_ued[par->index_ued[i][j][k]])/2.0;
			}
		}
	}
	
    
	return;
}


/* ******************************************************************* */
/*
 rearrange velocity u[] (j-2 & j+1 values for j)
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 
 u_cn is not defined at j = -1, ny+1
 */
/* ******************************************************************* */
void interpolate_y3_ued2cn__4(
                              double *u_ued,
                              double *u_cn,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				u_cn[par->index_cn[i][j][k]] \
                = (u_ued[par->index_ued[i][j-2][k]] + u_ued[par->index_ued[i][j+1][k]])/2.0;
			}
		}
	}
	
	//j = -1;
	//j = ny+1;
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			u_cn[par->index_cn[i][-1][k]] = 0.0;
			u_cn[par->index_cn[i][ny+1][k]] = 0.0;
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity u[] (index_u[][][]) (i & i+1 values for i)
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 u_ct is not defined at i = end+2
 */
/* ******************************************************************* */
void interpolate_x1_ued2ct__4(
                              double *u_ued,
                              double *u_ct,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start-3; i<=end+1; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_ct[par->index_ct[i][j][k]] \
                = (u_ued[par->index_ued[i][j][k]] + u_ued[par->index_ued[i+1][j][k]])/2.0;
			}
		}
	}
	
	i = end+2;
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			u_ct[par->index_ct[i][j][k]] = 0.0;
		}
	}
	
	return;
}



/* ******************************************************************* */
/*
 interpolate velocity u[] (index_u[][][]) (i-1 & i+2 values for i)
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 u_ct is not defined at i = start-3, end+1, end+2
 */
/* ******************************************************************* */
void interpolate_x3_ued2ct__4(
                              double *u_ued,
                              double *u_ct,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start-2; i<=end; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_ct[par->index_ct[i][j][k]] \
                = (u_ued[par->index_ued[i-1][j][k]] + u_ued[par->index_ued[i+2][j][k]])/2.0;
			}
		}
	}
	
	//i = start-3;
	//i = end+1;
	//i = end+2;
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			u_ct[par->index_ct[start-3][j][k]] = 0.0;
			u_ct[par->index_ct[end+1][j][k]] = 0.0;
			u_ct[par->index_ct[end+2][j][k]] = 0.0;
		}
	}
	
    
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity w[] (j & j-1 values for j)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_v[]: [start->end-1][1->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_y1_ct2v__4(
                            double *w_ct,
                            double *w_v,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_v[par->index_v[i][j][k]] \
                = (w_ct[par->index_ct[i][j-1][k]] + w_ct[par->index_ct[i][j][k]])/2.0;
			}
		}
	}
	
	
	return;
}




/* ******************************************************************* */
/*
 interpolate velocity w[] (j+1 & j-2 values for j)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_v[]: [start->end-1][1->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_y3_ct2v__4(
                            double *w_ct,
                            double *w_v,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_v[par->index_v[i][j][k]] \
                = (w_ct[par->index_ct[i][j-2][k]] + w_ct[par->index_ct[i][j+1][k]])/2.0;
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity w[] (i-1 & i values for i)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_u[]: [start->end-1][0->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_x1_ct2u__4(
                            double *w_ct,
                            double *w_u,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_u[par->index_u[i][j][k]] \
                = (w_ct[par->index_ct[i-1][j][k]] + w_ct[par->index_ct[i][j][k]])/2.0;
			}
		}
	}
	
	
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity w[] (i-2 & i+1 values for i)
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_u[]: [start->end-1][0->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_x3_ct2u__4(
                            double *w_ct,
                            double *w_u,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_u[par->index_u[i][j][k]] \
                = (w_ct[par->index_ct[i-2][j][k]] + w_ct[par->index_ct[i+1][j][k]])/2.0;
			}
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_v[]: [start->end-1][1->ny-1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_x1_cn2v__4(
                            double *v_cn,
                            double *v,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				v[par->index_v[i][j][k]] \
                = (v_cn[par->index_cn[i][j][k]] + v_cn[par->index_cn[i+1][j][k]])/2.0;
			}
		}
	}
	
	
	return;
}




/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_v[]: [start->end-1][1->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_x3_cn2v__4(
                            double *v_cn,
                            double *v,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				v[par->index_v[i][j][k]] \
                = (v_cn[par->index_cn[i-1][j][k]] + v_cn[par->index_cn[i+2][j][k]])/2.0;
			}
		}
	}
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 
 ved is not defined
 at j = -2, ny+2
 at i = start-3, start-2, end+1, end+2
 */
/* ******************************************************************* */
void interpolate_x1_cn2ved__4(
                              double *v_cn,
                              double *ved,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] \
				= (v_cn[par->index_cn[i][j][k]] + v_cn[par->index_cn[i+1][j][k]])/2.0;
			}
		}
	}
	
	
	for (i=start-3; i<=start-2; ++i){
		for (j=-2; j<=ny+2; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] = 0.0;
			}
		}
	}
	
	for (i=end+1; i<=end+2; ++i){
		for (j=-2; j<=ny+2; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] = 0.0;
			}
		}
	}
	
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -2;
			ved[par->index_ved[i][j][k]] = 0.0;
			j = ny+2;
			ved[par->index_ved[i][j][k]] = 0.0;
		}
	}
	
	return;
}




/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 
 ved is not defined
 at j = -2, ny+2
 at i = start-3, start-2, start-1
 at i = end, end+1, end+2
 */
/* ******************************************************************* */
void interpolate_x3_cn2ved__4(
                              double *v_cn,
                              double *ved,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] \
				= (v_cn[par->index_cn[i-1][j][k]] + v_cn[par->index_cn[i+2][j][k]])/2.0;
			}
		}
	}
	
	for (i=start-3; i<=start-1; ++i){
		for (j=-2; j<=ny+2; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] = 0.0;
			}
		}
	}
	
	for (i=end; i<=end+2; ++i){
		for (j=-2; j<=ny+2; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] = 0.0;
			}
		}
	}
	
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -2;
			ved[par->index_ved[i][j][k]] = 0.0;
			j = ny+2;
			ved[par->index_ved[i][j][k]] = 0.0;
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_u[]: [start->end-1][0->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_y1_cn2u__4(
                            double *u_cn,
                            double *u,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u[par->index_u[i][j][k]] \
                = (u_cn[par->index_cn[i][j][k]] + u_cn[par->index_cn[i][j+1][k]])/2.0;
			}
		}
	}
	
	
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_u[]: [start->end-1][0->ny-1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_y3_cn2u__4(
                            double *u_cn,
                            double *u,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u[par->index_u[i][j][k]] \
                = (u_cn[par->index_cn[i][j-1][k]] + u_cn[par->index_cn[i][j+2][k]])/2.0;
			}
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 i = start-3, start-2, end+2 will not be defined
 j = -2, ny+1 will not be defined
 */
/* ******************************************************************* */
void interpolate_y1_cn2ued__4(
                              double *u_cn,
                              double *u_ed,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				u_ed[par->index_ued[i][j][k]] \
				= (u_cn[par->index_cn[i][j][k]] + u_cn[par->index_cn[i][j+1][k]])/2.0;
			}
		}
	}
	
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			i = start-3;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			i = start-2;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			i = end+2;
			u_ed[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -2;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			j = ny+1;
			u_ed[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 i = start-3, start-2, end+2 will not be defined
 j = -1, -2, ny, ny+1 will not be defined
 */
/* ******************************************************************* */
void interpolate_y3_cn2ued__4(
                              double *u_cn,
                              double *u_ed,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_ed[par->index_ued[i][j][k]] \
				= (u_cn[par->index_cn[i][j-1][k]] + u_cn[par->index_cn[i][j+2][k]])/2.0;
			}
		}
	}
	
	
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			i = start-3;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			i = start-2;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			i = end+2;
			u_ed[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -1;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			j = -2;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			j = ny;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			j = ny+1;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity v[] (index_v[][][]) (j & j+1 values for j)
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_y1_ved2wp__4(
                              double *v_ved,
                              double *v_wp1,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				v_wp1[par->index_wp[i][j][k]] \
                = (v_ved[par->index_ved[i][j][k]] + v_ved[par->index_ved[i][j+1][k]])/2.0;
			}
		}
	}
	
	return;
}


/* ******************************************************************* */
/*
 interpolate velocity v[] (index_v[][][]) (j-1 & j+2 values for j)
 index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_y3_ved2wp__4(
                              double *v_ved,
                              double *v_wp,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				v_wp[par->index_wp[i][j][k]] \
                = (v_ved[par->index_ved[i][j-1][k]] + v_ved[par->index_ved[i][j+2][k]])/2.0;
			}
		}
	}
    
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity u[] (index_u[][][]) (i & i+1 values for i)
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_x1_ued2wp__4(
                              double *u_ued,
                              double *u_wp1,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_wp1[par->index_wp[i][j][k]] \
                = (u_ued[par->index_ued[i][j][k]] + u_ued[par->index_ued[i+1][j][k]])/2.0;
			}
		}
	}
    
	return;
}

/* ******************************************************************* */
/*
 interpolate velocity u[] (index_u[][][]) (i-1 & i+2 values for i)
 index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 
 */
/* ******************************************************************* */
void interpolate_x3_ued2wp__4(
                              double *u_ued,
                              double *u_wp3,
                              parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_wp3[par->index_wp[i][j][k]] \
                = (u_ued[par->index_ued[i-1][j][k]] + u_ued[par->index_ued[i+2][j][k]])/2.0;
			}
		}
	}
	
	
    
	return;
}

void interpolate_ct2cn__4(
						  double *ct,
						  double *cn,
						  parameters *par)
{
	int i;
	const int local_size_cn = par->local_size_cn;
	const int local_size_ved = par->local_size_ved;
	
	double *scratch1;
	double *scratch2;
	
	// interpolate ct onto v_edge ******************
	scratch1 = FFTW_MALLOC(local_size_ved);
	scratch2 = FFTW_MALLOC(local_size_ved);
	
	interpolate_y1_ct2ved__4(ct, scratch1, par);
	interpolate_y3_ct2ved__4(ct, scratch2, par);
	
	for (i=0; i<=local_size_ved-1; ++i){
		scratch1[i] = ( 9.0*scratch1[i] - 1.0*scratch2[i] )/8.0;
	}
	
	
	fftw_free(scratch2);
	
	
	
	// further interpolate onto corner ******************
	scratch2 = FFTW_MALLOC(local_size_cn);
	
	interpolate_x1_ved2cn__4(scratch1, cn, par);
	interpolate_x3_ved2cn__4(scratch1, scratch2, par);
	
	for (i=0; i<=local_size_cn-1; ++i){
		cn[i] = ( 9.0*cn[i] - 1.0*scratch2[i] )/8.0;
	}
	
	fftw_free(scratch1);
	fftw_free(scratch2);
	
	
	return;
}



/* ******************************************************************* */
/*
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 */
/* ******************************************************************* */

void interpolate_cn2wp__4(
						  double *cn,
						  double *wp,
						  parameters *par)
{
	int i;
	const int local_size_wp = par->local_size_wp;
	const int local_size_ued = par->local_size_ued;
	
	double *scratch1;
	double *scratch2;
	
	// interpolate dudy onto u_edge ******************
	// index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]
	// index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
	// i = start-3, start-2, end+2 will not be defined
	// j = -1, -2, ny, ny+1 will not be defined
	
	scratch1 = FFTW_MALLOC(local_size_ued);
	scratch2 = FFTW_MALLOC(local_size_ued);
	
	interpolate_y1_cn2ued__4(cn, scratch1, par);
	interpolate_y3_cn2ued__4(cn, scratch2, par);
	
	for (i=0; i<=local_size_ued-1; ++i){
		scratch1[i] = ( 9.0*scratch1[i] - 1.0*scratch2[i] )/8.0;
	}
	
	fftw_free(scratch2);
	
	
	
	// further interpolated onto center ******************
	// For scratch1,
	// index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]
	// i = start-3, start-2, end+2 will not be defined
	// j = -1, -2, ny, ny+1 will not be defined
	// But enough to define interpolated wp
	
	scratch2 = FFTW_MALLOC(local_size_wp);
	
	interpolate_x1_ued2wp__4(scratch1, wp, par);
	interpolate_x3_ued2wp__4(scratch1, scratch2, par);
	
	for (i=0; i<=local_size_wp-1; ++i){
		wp[i] = ( 9.0*wp[i] - 1.0*scratch2[i] )/8.0;
	}
	
	fftw_free(scratch1);
	fftw_free(scratch2);
	
	
	return;
}

/* ******************************************************************* */
/*
 rearrange velocity u[] (index_u[][][])
 
 index_ued[] is defined at the edges of each cell similar to u[], but boundary
 values are included.
 index_ued[]: [-width->nx-1+width][-width_uw->ny-1+width_uw][0->nz-1]
 
 flag = 0 (homogeneous boundary)
 flag = 1 (inhomogeneous boudary)
 */
/* ******************************************************************* */
void bigger_array_u_ued__4(
                           double flag,
                           variables *var,
                           neighbors *shared,
                           double *u_ed,
                           parameters *par)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int width = par->width_share;
	const int width_uw = 2;
	
	double *u_r;
	double *u_l;
    
	int ***index_share_r;
	int ***index_share_l;
	
	index_share_l = i3tensor(start-width, start-1, -width_uw, ny-1+width_uw, 0, nz-1);
	index_share_r = i3tensor(end, end-1+width, -width_uw, ny-1+width_uw, 0, nz-1);
	
	for (i=start-width; i<=start-1; ++i){
		for (j=-width_uw; j<=ny-1+width_uw; ++j){
			for (k=0; k<=nz-1; ++k){
				index_share_l[i][j][k] = (i-start+width)*((ny-1+2*width_uw)*nz + (nz-1) + 1) + (j+width_uw)*nz + k;
			}
		}
	}
	
	for (i=end; i<=end-1+width; ++i){
		for (j=-width_uw; j<=ny-1+width_uw; ++j){
			for (k=0; k<=nz-1; ++k){
				index_share_r[i][j][k] = (i-end)*((ny-1+2*width_uw)*nz + (nz-1) + 1) + (j+width_uw)*nz + k;
			}
		}
	}
	
    
	u_r = dvector(0, width*(ny+2*width_uw)*nz-1);
	u_l = dvector(0, width*(ny+2*width_uw)*nz-1);
	
	for (i=0; i<=width*(ny+2*width_uw)*nz-1; ++i){
		u_r[i] = 0.0;
		u_l[i] = 0.0;
	}
    
	index = 0;
	for (i=start-width; i<=start-1; ++i){
		// bottom ghost cell
		for (k=0; k<=nz-1; ++k){
			u_l[index_share_l[i][-2][k]] = 8.0*flag*shared->u_in_l[index] \
            - 9.0*shared->u_in_l[index+nz] \
            + 2.0*shared->u_in_l[index+2*nz];
            
			u_l[index_share_l[i][-1][k]] = 8.0/3.0*flag*shared->u_in_l[index] \
            - 2.0*shared->u_in_l[index+nz] \
            + 1.0/3.0*shared->u_in_l[index+2*nz];
			index = index + 1;
		}
		// interior
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_l[index_share_l[i][j][k]] = shared->u_in_l[index];
				index = index + 1;
			}
		}
		// top ghost cell
		for (k=0; k<=nz-1; ++k){
			u_l[index_share_l[i][ny][k]] = 8.0/3.0*flag*shared->u_in_l[index] \
            - 2.0*shared->u_in_l[index-nz] \
            + 1.0/3.0*shared->u_in_l[index-2*nz];
			
			u_l[index_share_l[i][ny+1][k]] = 8.0*flag*shared->u_in_l[index] \
            - 9.0*shared->u_in_l[index-nz] \
            + 2.0*shared->u_in_l[index-2*nz];
			index = index + 1;
		}
	}
    
	
	index = 0;
	for (i=end; i<=end-1+width; ++i){
		// bottom ghost cell
		for (k=0; k<=nz-1; ++k){
			u_r[index_share_r[i][-2][k]] = 8.0*flag*shared->u_in_r[index] \
            - 9.0*shared->u_in_r[index+nz] \
            + 2.0*shared->u_in_r[index+2*nz];
            
			u_r[index_share_r[i][-1][k]] = 8.0/3.0*flag*shared->u_in_r[index] \
            - 2.0*shared->u_in_r[index+nz] \
            + 1.0/3.0*shared->u_in_r[index+2*nz];
			index = index + 1;
		}
		// interior
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_r[index_share_r[i][j][k]] = shared->u_in_r[index];
				index = index + 1;
			}
		}
		// top ghost cell
		for (k=0; k<=nz-1; ++k){
			u_r[index_share_r[i][ny][k]] = 8.0/3.0*flag*shared->u_in_r[index] \
            - 2.0*shared->u_in_r[index-nz] \
            + 1.0/3.0*shared->u_in_r[index-2*nz];
			
			u_r[index_share_r[i][ny+1][k]] = 8.0*flag*shared->u_in_r[index] \
            - 9.0*shared->u_in_r[index-nz] \
            + 2.0*shared->u_in_r[index-2*nz];
			index = index + 1;
		}
	}
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// index_ued[]: [-1->nx][-1->ny][0->nz-1]
	//  [-width->nx-1+width][-width_uw->ny-1+width_uw][0->nz-1]
	// define ued at the edges
	for (i=start-width; i<=end-1+width; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				// left boundary
				if (i < start){
					u_ed[par->index_ued[i][j][k]] = u_l[index_share_l[i][j][k]];
                    // right boundary
				}else if (i > end-1){
					u_ed[par->index_ued[i][j][k]] = u_r[index_share_r[i][j][k]];
                    // interior
				}else{
					u_ed[par->index_ued[i][j][k]] = var->u[par->index_u[i][j][k]];
				}
			}
		}
	}
	
	// bottom boundary
	for (i=start-width; i<=end-1+width; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i < start){
				u_ed[par->index_ued[i][-2][k]] = u_l[index_share_l[i][-2][k]];
				u_ed[par->index_ued[i][-1][k]] = u_l[index_share_l[i][-1][k]];
                // right boundary
			}else if (i > end-1){
				u_ed[par->index_ued[i][-2][k]] = u_r[index_share_r[i][-2][k]];
				u_ed[par->index_ued[i][-1][k]] = u_r[index_share_r[i][-1][k]];
                // interior
			}else{
				u_ed[par->index_ued[i][-1][k]] = 8.0/3.0*flag*var->bc_u_b[par->index_tb[i][k]] \
                - 2.0*var->u[par->index_u[i][0][k]]	\
                +1.0/3.0*var->u[par->index_u[i][1][k]];
				u_ed[par->index_ued[i][-2][k]] = 8.0*flag*var->bc_u_b[par->index_tb[i][k]] \
                - 9.0*var->u[par->index_u[i][0][k]]	\
                + 2.0*var->u[par->index_u[i][1][k]];
			}
		}
	}
    
	
	// top boundary
	for (i=start-width; i<=end-1+width; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i < start){
				u_ed[par->index_ued[i][ny][k]] = u_l[index_share_l[i][ny][k]];
				u_ed[par->index_ued[i][ny+1][k]] = u_l[index_share_l[i][ny+1][k]];
                // right boundary
			}else if (i > end-1){
				u_ed[par->index_ued[i][ny][k]] = u_r[index_share_r[i][ny][k]];
				u_ed[par->index_ued[i][ny+1][k]] = u_r[index_share_r[i][ny+1][k]];
                // interior
			}else{
				u_ed[par->index_ued[i][ny][k]] = 8.0/3.0*flag*var->bc_u_t[par->index_tb[i][k]] \
                - 2.0*var->u[par->index_u[i][ny-1][k]]	\
                +1.0/3.0*var->u[par->index_u[i][ny-2][k]];
				u_ed[par->index_ued[i][ny+1][k]] = 8.0*flag*var->bc_u_t[par->index_tb[i][k]] \
                - 9.0*var->u[par->index_u[i][ny-1][k]]	\
                + 2.0*var->u[par->index_u[i][ny-2][k]];
			}
		}
	}
    
	
	free_dvector(u_r, 0, width*(ny+2*width_uw)*nz-1);
	free_dvector(u_l, 0, width*(ny+2*width_uw)*nz-1);
	
	free_i3tensor(index_share_l, start-width, start-1, -width_uw, ny-1+width_uw, 0, nz-1);
	free_i3tensor(index_share_r, end, end-1+width, -width_uw, ny-1+width_uw, 0, nz-1);
    
	return;
	
}

/* ******************************************************************* */
/*
 rearrange velocity v[] (index_v[][][])
 
 index_ved[] is defined at the top edge of each cell similar to v[], but boundary
 values are included. index_ved[]:
 (start-width, end-1+width, -width_v, ny-1+width_v, 0, nz-1);
 
 y = -width_v (-2) and y = ny-1+width_v (ny+2) are to be defined later
 for calculating non linear terms
 those values are defined b/w i = -2 to i = end
 by function bigger_array_v_ved__4_addition()
 
 
 flag = 0 (homogeneous boundary)
 flag = 1 (inhomogeneous boudary)
 */
/* ******************************************************************* */
void bigger_array_v_ved__4(
                           double flag,
                           variables *var,
                           neighbors *shared,
                           double *v_ed,
                           parameters *par)
{
	int i, j, k, index;
	int my_rank; my_rank = par->my_rank;
	int nx; nx = par->nx;
	int ny; ny = par->ny;
	int nz; nz = par->nz;
	int start = par->local_nx_start;
	int end = par->local_nx_start + par->local_nx;
	int width; width = par->width_share;
	
	double *v_r;
	double *v_l;
	
	int ***index_share_r;
	int ***index_share_l;
	
	index_share_l = i3tensor(start-width, start-1, -1, ny+1, 0, nz-1);
	index_share_r = i3tensor(end, end-1+width, -1, ny+1, 0, nz-1);
	
	for (i=start-width; i<=start-1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				index_share_l[i][j][k] = (i-start+width)*((ny+2)*nz + (nz-1) + 1) + (j+1)*nz + k;
			}
		}
	}
	
	for (i=end; i<=end-1+width; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				index_share_r[i][j][k] = (i-end)*((ny+2)*nz + (nz-1) + 1) + (j+1)*nz + k;
			}
		}
	}
	
	v_r = dvector(0, width*(ny+3)*nz-1);
	v_l = dvector(0, width*(ny+3)*nz-1);
	
	for (i=0; i<=width*(ny+3)*nz-1; ++i){
		v_r[i] = 0.0;
		v_l[i] = 0.0;
	}
	
	index = 0;
	for (i=start-width; i<=start-1; ++i){
		// bottom ghost cell
		for (k=0; k<=nz-1; ++k){
			v_l[index_share_l[i][-1][k]] = 2.0*flag*shared->v_in_l[index] - shared->v_in_l[index+nz];
            index = index + 1;
		}
		index = index - nz;
		// interior
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				if (flag == 0.0 && (j == 0 || j == ny)){
					v_l[index_share_l[i][j][k]] = flag*shared->v_in_l[index];
				}else{
					v_l[index_share_l[i][j][k]] = shared->v_in_l[index];
				}
				index = index + 1;
			}
		}
		index = index - nz;
		// top ghost cell
		for (k=0; k<=nz-1; ++k){
			v_l[index_share_l[i][ny+1][k]] = 2.0*flag*shared->v_in_l[index] - shared->v_in_l[index-nz];
			index = index + 1;
		}
	}
	
	index = 0;
	for (i=end; i<=end-1+width; ++i){
		// bottom ghost cell
		for (k=0; k<=nz-1; ++k){
			v_r[index_share_r[i][-1][k]] \
            = 2.0*flag*shared->v_in_r[index] \
            - shared->v_in_r[index+nz];
			index = index + 1;
		}
		index = index - nz;
		// interior
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				if (flag == 0.0 && (j == 0 || j == ny)){
					v_r[index_share_r[i][j][k]] = flag*shared->v_in_r[index];
				}else{
					v_r[index_share_r[i][j][k]] = shared->v_in_r[index];
				}
				index = index + 1;
			}
		}
		index = index - nz;
		// top ghost cell
		for (k=0; k<=nz-1; ++k){
			v_r[index_share_r[i][ny+1][k]] \
            = 2.0*flag*shared->v_in_r[index] \
            - shared->v_in_r[index-nz];
			index = index + 1;
		}
	}
	
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// index_ved[]: (start-width, end-1+width, 1-width_v, ny-1+width_v, 0, nz-1)
	// define ued at the edges
	for (i=start-width; i<=end-1+width; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				// left boundary
				if (i < start){
					v_ed[par->index_ved[i][j][k]] = v_l[index_share_l[i][j][k]];
                    // right boundary
				}else if (i > end-1){
					v_ed[par->index_ved[i][j][k]] = v_r[index_share_r[i][j][k]];
                    // interior
				}else{
					v_ed[par->index_ved[i][j][k]] = var->v[par->index_v[i][j][k]];
				}
			}
		}
	}
	
	for (i=start-width; i<=end-1+width; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i < start){
				v_ed[par->index_ved[i][0][k]] = v_l[index_share_l[i][0][k]];
				v_ed[par->index_ved[i][-1][k]] = v_l[index_share_l[i][-1][k]];
                // right boundary
			}else if (i > end-1){
				v_ed[par->index_ved[i][0][k]] = v_r[index_share_r[i][0][k]];
				v_ed[par->index_ved[i][-1][k]] = v_r[index_share_r[i][-1][k]];
                // interior
			}else{
				v_ed[par->index_ved[i][0][k]]
                = flag*var->bc_v_b[par->index_tb[i][k]];
				v_ed[par->index_ved[i][-1][k]]
                = 2.0*flag*var->bc_v_b[par->index_tb[i][k]]
                - var->v[par->index_v[i][1][k]];
                
			}
		}
	}
	
	for (i=start-width; i<=end-1+width; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i < start){
				v_ed[par->index_ved[i][ny][k]] = v_l[index_share_l[i][ny][k]];
				v_ed[par->index_ved[i][ny+1][k]] = v_l[index_share_l[i][ny+1][k]];
                // right boundary
			}else if (i > end-1){
				v_ed[par->index_ved[i][ny][k]] = v_r[index_share_r[i][ny][k]];
				v_ed[par->index_ved[i][ny+1][k]] = v_r[index_share_r[i][ny+1][k]];
                // interior
			}else{
				v_ed[par->index_ved[i][ny][k]] = flag*var->bc_v_t[par->index_tb[i][k]];
				v_ed[par->index_ved[i][ny+1][k]] \
                = 2.0*flag*var->bc_v_t[par->index_tb[i][k]] \
                - var->v[par->index_v[i][ny-1][k]];
            }
		}
	}
	
    
	free_dvector(v_r, 0, width*(ny+3)*nz-1);
	free_dvector(v_l, 0, width*(ny+3)*nz-1);
	
    
	free_i3tensor(index_share_l, start-width, start-1, -1, ny+1, 0, nz-1);
	free_i3tensor(index_share_r, end, end-1+width, -1, ny+1, 0, nz-1);
	
	return;
}

/* ******************************************************************* */
/*
 y = -width_v (-2) and y = ny-1+width_v (ny+2) are to be defined.
 these values are defined b/w i = start-2 to i = end
 by function bigger_array_v_ved__4_addition (here)
 */
/* ******************************************************************* */
void bigger_array_v_ved__4_addition(
                                    double *u_ed,
                                    double *v_ed,
                                    double *w_ct,
                                    parameters *par)
{
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
    const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	double *w_t, *w_b;
	double *dwdz_t, *dwdz_b;
	double *dudx_t, *dudx_b;
	
	dwdz_t = FFTW_MALLOC((local_nx+3)*nz);
	dwdz_b = FFTW_MALLOC((local_nx+3)*nz);
	dudx_t = FFTW_MALLOC((local_nx+3)*nz);
	dudx_b = FFTW_MALLOC((local_nx+3)*nz);
	
	w_t = FFTW_MALLOC((local_nx+3)*nz);
	w_b = FFTW_MALLOC((local_nx+3)*nz);
    
	for (i=start-2; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			w_b[(i-start+2)*nz + k] = w_ct[par->index_ct[i][-1][k]];
			w_t[(i-start+2)*nz + k] = w_ct[par->index_ct[i][ny][k]];
            
		}
	}
	
	for (i=start-2; i<=end; ++i){
		k = 0;
		dwdz_b[(i-start+2)*nz + k] = 0.0;
		dwdz_t[(i-start+2)*nz + k] = 0.0;
        
		if (nz%2 == 0){
			k = nz/2;
			dwdz_b[(i-start+2)*nz + k] = 0.0;
			dwdz_t[(i-start+2)*nz + k] = 0.0;
		}
        
		for (k=1; k<=(nz-1)/2; ++k){
			dwdz_b[(i-start+2)*nz + (nz-k)] \
            = w_b[(i-start+2)*nz + k]*par->kz[k];
			dwdz_b[(i-start+2)*nz + k] \
            = - w_b[(i-start+2)*nz + (nz-k)]*par->kz[nz-k];
			
			dwdz_t[(i-start+2)*nz + (nz-k)] \
            = w_t[(i-start+2)*nz + k]*par->kz[k];
			dwdz_t[(i-start+2)*nz + k] \
            = - w_t[(i-start+2)*nz + (nz-k)]*par->kz[nz-k];
		}
	}
	
	fftw_free(w_b);
	fftw_free(w_t);
    
	
	
	for (i=start-2; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -1;
			dudx_b[(i-start+2)*nz + k] \
			= (u_ed[par->index_ued[i-1][j][k]] - 27.0*u_ed[par->index_ued[i][j][k]] \
               + 27.0*u_ed[par->index_ued[i+1][j][k]] - u_ed[par->index_ued[i+2][j][k]])/(24.0*par->dx);
			
			j = ny;
			dudx_t[(i-start+2)*nz + k] \
			= (u_ed[par->index_ued[i-1][j][k]] - 27.0*u_ed[par->index_ued[i][j][k]] \
               + 27.0*u_ed[par->index_ued[i+1][j][k]] - u_ed[par->index_ued[i+2][j][k]])/(24.0*par->dx);
		}
	}
    
    
	for (i=start-2; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			v_ed[par->index_ved[i][-2][k]] \
            = 4.0*( v_ed[par->index_ved[i][1][k]] + v_ed[par->index_ved[i][-1][k]] ) \
            - v_ed[par->index_ved[i][2][k]] \
            - 6.0*v_ed[par->index_ved[i][0][k]];
            v_ed[par->index_ved[i][ny+2][k]]
            = 4.0*( v_ed[par->index_ved[i][ny-1][k]] + v_ed[par->index_ved[i][ny+1][k]] )
            - v_ed[par->index_ved[i][ny-2][k]]
            - 6.0*v_ed[par->index_ved[i][ny][k]];
		}
	}
    
	
	fftw_free(dwdz_b);
	fftw_free(dwdz_t);
	fftw_free(dudx_b);
	fftw_free(dudx_t);
	
	
	//j = -2;
	//j = ny+2:
	i = start-3;
	for (k=0; k<=nz-1; ++k){
		v_ed[par->index_ved[i][-2][k]] = 0.0;
		v_ed[par->index_ved[i][ny+2][k]] = 0.0;
	}
	
	i = end+2;
	for (k=0; k<=nz-1; ++k){
		v_ed[par->index_ved[i][-2][k]] = 0.0;
		v_ed[par->index_ved[i][ny+2][k]] = 0.0;
	}
	
	i = end+1;
	for (k=0; k<=nz-1; ++k){
		v_ed[par->index_ved[i][-2][k]] = 0.0;
		v_ed[par->index_ved[i][ny+2][k]] = 0.0;
	}
	
	return;
}



/* ******************************************************************* */
/*
 rearrange velocity w[] (index_wp[][][])
 
 index_ct[] is defined as w[], but boundary values are included.
 
 index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]
 
 flag = 0 (homogeneous boundary)
 flag = 1 (inhomogeneous boudary)
 */
/* ******************************************************************* */
void bigger_array_wp_ct__4(
                           double flag,
                           variables *var,
                           neighbors *shared,
                           double *w_ct,
                           parameters *par)
{
	
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int width = par->width_share;
	const int width_uw = 2;
	
	double *w_r;
	double *w_l;
    
	int ***index_share_r;
	int ***index_share_l;
	
	index_share_l = i3tensor(start-width, start-1, -width_uw, ny-1+width_uw, 0, nz-1);
	index_share_r = i3tensor(end, end-1+width, -width_uw, ny-1+width_uw, 0, nz-1);
	
	for (i=start-width; i<=start-1; ++i){
		for (j=-width_uw; j<=ny-1+width_uw; ++j){
			for (k=0; k<=nz-1; ++k){
				index_share_l[i][j][k] = (i-start+width)*((ny-1+2*width_uw)*nz + (nz-1) + 1) + (j+width_uw)*nz + k;
			}
		}
	}
	
	for (i=end; i<=end-1+width; ++i){
		for (j=-width_uw; j<=ny-1+width_uw; ++j){
			for (k=0; k<=nz-1; ++k){
				index_share_r[i][j][k] = (i-end)*((ny-1+2*width_uw)*nz + (nz-1) + 1) + (j+width_uw)*nz + k;
			}
		}
	}
    
	
	w_r = dvector(0, width*(ny+2*width_uw)*nz-1);
	w_l = dvector(0, width*(ny+2*width_uw)*nz-1);
	
	for (i=0; i<=width*(ny+2*width_uw)*nz-1; ++i){
		w_r[i] = 0.0;
		w_l[i] = 0.0;
	}
    
	index = 0;
	for (i=start-width; i<=start-1; ++i){
		// bottom ghost cell
		for (k=0; k<=nz-1; ++k){
			w_l[index_share_l[i][-2][k]] = 8.0*flag*shared->w_in_l[index] \
            - 9.0*shared->w_in_l[index+nz] \
            + 2.0*shared->w_in_l[index+2*nz];
			
			w_l[index_share_l[i][-1][k]] = 8.0/3.0*flag*shared->w_in_l[index] \
            - 2.0*shared->w_in_l[index+nz] \
            + 1.0/3.0*shared->w_in_l[index+2*nz];
			
            
			index += 1;
		}
		// interior
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_l[index_share_l[i][j][k]] = shared->w_in_l[index];
				index = index + 1;
			}
		}
		// top ghost cell
		for (k=0; k<=nz-1; ++k){
			w_l[index_share_l[i][ny][k]] = 8.0/3.0*flag*shared->w_in_l[index] \
            - 2.0*shared->w_in_l[index-nz] \
            + 1.0/3.0*shared->w_in_l[index-2*nz];
			
			w_l[index_share_l[i][ny+1][k]] = 8.0*flag*shared->w_in_l[index] \
            - 9.0*shared->w_in_l[index-nz] \
            + 2.0*shared->w_in_l[index-2*nz];
			index = index + 1;
		}
	}
	
    
	index = 0;
	for (i=end; i<=end-1+width; ++i){
		// bottom ghost cell
		for (k=0; k<=nz-1; ++k){
			w_r[index_share_r[i][-2][k]] = 8.0*flag*shared->w_in_r[index] \
            - 9.0*shared->w_in_r[index+nz] \
            + 2.0*shared->w_in_r[index+2*nz];
			
			w_r[index_share_r[i][-1][k]] = 8.0/3.0*flag*shared->w_in_r[index] \
            - 2.0*shared->w_in_r[index+nz] \
            + 1.0/3.0*shared->w_in_r[index+2*nz];
            
			index += 1;
		}
		// interior
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_r[index_share_r[i][j][k]] = shared->w_in_r[index];
				index = index + 1;
			}
		}
		// top ghost cell
		for (k=0; k<=nz-1; ++k){
			w_r[index_share_r[i][ny][k]] = 8.0/3.0*flag*shared->w_in_r[index] \
            - 2.0*shared->w_in_r[index-nz] \
            + 1.0/3.0*shared->w_in_r[index-2*nz];
			
			w_r[index_share_r[i][ny+1][k]] = 8.0*flag*shared->w_in_r[index] \
            - 9.0*shared->w_in_r[index-nz] \
            + 2.0*shared->w_in_r[index-2*nz];
			index = index + 1;
		}
	}
    
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// index_ct[]: [-width->nx-1+width][-width_uw->ny-1+width_uw][0->nz-1]
	// define ued at the edges
	for (i=start-width; i<=end-1+width; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				// left boundary
				if (i < start){
					w_ct[par->index_ct[i][j][k]] = w_l[index_share_l[i][j][k]];
                    // right boundary
				}else if (i > end-1){
					w_ct[par->index_ct[i][j][k]] = w_r[index_share_r[i][j][k]];
                    // interior
				}else{
					w_ct[par->index_ct[i][j][k]] = var->w[par->index_wp[i][j][k]];
				}
			}
		}
	}
	
	// bottom boundary
	for (i=start-width; i<=end-1+width; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i < start){
				w_ct[par->index_ct[i][-2][k]] = w_l[index_share_l[i][-2][k]];
				w_ct[par->index_ct[i][-1][k]] = w_l[index_share_l[i][-1][k]];
                // right boundary
			}else if (i > end-1){
				w_ct[par->index_ct[i][-2][k]] = w_r[index_share_r[i][-2][k]];
				w_ct[par->index_ct[i][-1][k]] = w_r[index_share_r[i][-1][k]];
                // interior
			}else{
				w_ct[par->index_ct[i][-1][k]] = 8.0/3.0*flag*var->bc_w_b[par->index_tb[i][k]] \
                - 2.0*var->w[par->index_wp[i][0][k]]	\
                + 1.0/3.0*var->w[par->index_wp[i][1][k]];
				w_ct[par->index_ct[i][-2][k]] = 8.0*flag*var->bc_w_b[par->index_tb[i][k]] \
                - 9.0*var->w[par->index_wp[i][0][k]]	\
                + 2.0*var->w[par->index_wp[i][1][k]];
			}
		}
	}
    
    
	
	
	
	// top boundary
	for (i=start-width; i<=end-1+width; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i < start){
				w_ct[par->index_ct[i][ny][k]] = w_l[index_share_l[i][ny][k]];
				w_ct[par->index_ct[i][ny+1][k]] = w_l[index_share_l[i][ny+1][k]];
                // right boundary
			}else if (i > end-1){
				w_ct[par->index_ct[i][ny][k]] = w_r[index_share_r[i][ny][k]];
				w_ct[par->index_ct[i][ny+1][k]] = w_r[index_share_r[i][ny+1][k]];
                // interior
			}else{
				w_ct[par->index_ct[i][ny][k]] = 8.0/3.0*flag*var->bc_w_t[par->index_tb[i][k]] \
                - 2.0*var->w[par->index_wp[i][ny-1][k]]	\
                + 1.0/3.0*var->w[par->index_wp[i][ny-2][k]];
				w_ct[par->index_ct[i][ny+1][k]] = 8.0*flag*var->bc_w_t[par->index_tb[i][k]] \
                - 9.0*var->w[par->index_wp[i][ny-1][k]]	\
                + 2.0*var->w[par->index_wp[i][ny-2][k]];
			}
		}
	}
	
	
	free_dvector(w_r, 0, width*(ny+2*width_uw)*nz-1);
	free_dvector(w_l, 0, width*(ny+2*width_uw)*nz-1);
	
	free_i3tensor(index_share_l, start-width, start-1, -width_uw, ny-1+width_uw, 0, nz-1);
	free_i3tensor(index_share_r, end, end-1+width, -width_uw, ny-1+width_uw, 0, nz-1);
    
	return;
}

