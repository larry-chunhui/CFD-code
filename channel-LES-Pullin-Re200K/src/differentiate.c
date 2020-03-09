/*
 *  differentiate.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Differentiations with second and fourth schemes.
 */

#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "differentiate.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

/* ************************************************************************ */
/* ************************************************************************ */
void derivative_x1_ct2u(
                        double *uu1,
                        double *scratch,
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
				scratch[par->index_u[i][j][k]] \
				= (uu1[par->index_ct[i][j][k]] - uu1[par->index_ct[i-1][j][k]])/(par->dx);
			}
		}
	}
    
	return;
}



/* ************************************************************************ */
/* ************************************************************************ */
void derivative_x3_ct2u(
                        double *uu3,
                        double *scratch,
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
				scratch[par->index_u[i][j][k]] \
				= (uu3[par->index_ct[i+1][j][k]] - uu3[par->index_ct[i-2][j][k]])/(3.0*par->dx);
			}
		}
	}
	
	return;
}



/* ************************************************************************ */
/* ************************************************************************ */
void derivative_x1_cn2v(
                        double *uv1,
                        double *scratch,
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
				scratch[par->index_v[i][j][k]] \
                = (uv1[par->index_cn[i+1][j][k]] - uv1[par->index_cn[i][j][k]])/(par->dx);
			}
		}
	}
	
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_x3_cn2v(
                        double *uv3,
                        double *scratch,
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
				scratch[par->index_v[i][j][k]] \
                = (uv3[par->index_cn[i+2][j][k]] - uv3[par->index_cn[i-1][j][k]])/(3.0*par->dx);
			}
		}
	}
	
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_x1_ued2wp(
                          double *uw1,
                          double *scratch,
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
				scratch[par->index_wp[i][j][k]] \
                = (uw1[par->index_ued[i+1][j][k]] - uw1[par->index_ued[i][j][k]])/(par->dx);
			}
		}
	}
	
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_x3_ued2wp(
                          double *uw3,
                          double *scratch,
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
				scratch[par->index_wp[i][j][k]] \
                = (uw3[par->index_ued[i+2][j][k]] - uw3[par->index_ued[i-1][j][k]])/(3.0*par->dx);
			}
		}
	}
	
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_y1_cn2u(
                        double *uv1,
                        double *scratch,
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
				scratch[par->index_u[i][j][k]] \
                = (uv1[par->index_cn[i][j+1][k]] - uv1[par->index_cn[i][j][k]])/(par->dy);
			}
		}
	}
	
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_y3_cn2u(
                        double *uv3,
                        double *scratch,
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
				scratch[par->index_u[i][j][k]] \
                = (uv3[par->index_cn[i][j+2][k]] - uv3[par->index_cn[i][j-1][k]])/(3.0*par->dy);
			}
		}
	}
	
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_y1_ct2v(
                        double *vv1,
                        double *scratch,
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
				scratch[par->index_v[i][j][k]] \
                =(vv1[par->index_ct[i][j][k]] - vv1[par->index_ct[i][j-1][k]])/(par->dy);
			}
		}
	}
    
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_y3_ct2v(
                        double *vv3,
                        double *scratch,
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
				scratch[par->index_v[i][j][k]] \
                =(vv3[par->index_ct[i][j+1][k]] - vv3[par->index_ct[i][j-2][k]])/(3.0*par->dy);
			}
		}
	}
    
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_y1_ved2wp(
                          double *vw1,
                          double *scratch,
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
				scratch[par->index_wp[i][j][k]] \
                = (vw1[par->index_ved[i][j+1][k]] - vw1[par->index_ved[i][j][k]])/(par->dy);
			}
		}
	}
    
	return;
}


/* ************************************************************************ */
/* ************************************************************************ */
void derivative_y3_ved2wp(
                          double *vw3,
                          double *scratch,
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
				scratch[par->index_wp[i][j][k]] \
                = (vw3[par->index_ved[i][j+2][k]] - vw3[par->index_ved[i][j-1][k]])/(3.0*par->dy);
			}
		}
	}
    
    
	return;
}


/* ************************************************************************ */
// d(uw)/dz defined at the same position is velocity u[] (edge)
// index_u[][][] is used as a location indicator
// while index_ued[][][] is for uw[]
/* ************************************************************************ */
void derivative_z_ued2u(
                        double *uw,
                        double *scratch,
                        parameters *par)
{
    
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			k = 0;
			scratch[par->index_u[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				scratch[par->index_u[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				scratch[par->index_u[i][j][nz-k]] = uw[par->index_ued[i][j][k]]*par->kz[k];
				scratch[par->index_u[i][j][k]] =  - uw[par->index_ued[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
    
	return;
}


/* ************************************************************************ */
// d(vw)/dz defined at the same position is velocity v[] (edge)
// index_v[][][] is used as a location indicator
// while index_ved[][][] is for vw[]
/* ************************************************************************ */
void derivative_z_ved2v(
                        double *vw,
                        double *scratch,
                        parameters *par)
{
    
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
    
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			k = 0;
			scratch[par->index_v[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				scratch[par->index_v[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				scratch[par->index_v[i][j][nz-k]] = vw[par->index_ved[i][j][k]]*par->kz[k];
				scratch[par->index_v[i][j][k]] =  - vw[par->index_ved[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
    
	return;
}




/* ************************************************************************ */
// d(ww)/dz defined at the same position is velocity w[] (center)
// index_wp[][][] is used as a location indicator
// same for ww[]
/* ************************************************************************ */
void derivative_z_ct2wp(
                        double *ww,
                        double *scratch,
                        parameters *par)
{
    
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			k = 0;
			scratch[par->index_wp[i][j][k]] = 0.0;
			
			if (nz%2 == 0){
				k = nz/2;
				scratch[par->index_wp[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				scratch[par->index_wp[i][j][nz-k]] = ww[par->index_ct[i][j][k]]*par->kz[k];
				scratch[par->index_wp[i][j][k]] =  - ww[par->index_ct[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
    
	return;
}


/* ******************************************************************* */
/* take a derivative in x (stream direction) for u (index_ued)		   */
/* derivatives are defined at center, index_ct						   */
/* index_ued[]: [start-1->end][-1->ny][0->nz-1]					   */
/* index_ct[]: [start-1->end][-1->ny][0->nz-1]					   */
/* u_ct is not defined at i = end									   */
/* ******************************************************************* */
void derivative_x1_ued2ct__2(
							 double *u,
							 double *dudx1,
							 parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const double dx = par->dx;
	
	//interior
	for (i=start-1; i<=end-1; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dudx1[par->index_ct[i][j][k]] \
				= (u[par->index_ued[i+1][j][k]] - u[par->index_ued[i][j][k]])/dx;
			}
		}
	}
	
	i = end;
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			dudx1[par->index_ct[i][j][k]] = 0.0;
		}
	}
	
	
	return;
}






/* ******************************************************************* */
/* take a derivative in x (stream direction) for u (index_ued)		   */
/* derivatives are defined at center, index_ct						   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* u_ct is not defined at i = end+2									   */
/* ******************************************************************* */
void derivative_x1_ued2ct__4(
                             double *u,
                             double *dudx1,
                             parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double dx = par->dx;
	
	//interior
	for (i=start-3; i<=end+1; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dudx1[par->index_ct[i][j][k]] \
                = (u[par->index_ued[i+1][j][k]] - u[par->index_ued[i][j][k]])/dx;
			}
		}
	}
	
	i = end+2;
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			dudx1[par->index_ct[i][j][k]] = 0.0;
		}
	}
    
    
	return;
}


/* ******************************************************************* */
/* take a derivative in x (stream direction) for u (index_ued)		   */
/* derivatives are defined at center, index_ct						   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* u_ct is not defined at i = start-3, end+2, end+1					   */
/* ******************************************************************* */
void derivative_x3_ued2ct__4(
                             double *u,
                             double *dudx,
                             parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double dx = par->dx;
	
	for (i=start-2; i<=end; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dudx[par->index_ct[i][j][k]] \
                = (u[par->index_ued[i+2][j][k]] - u[par->index_ued[i-1][j][k]])/(3.0*dx);
			}
		}
	}
	
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			dudx[par->index_ct[start-3][j][k]] = 0.0;
			dudx[par->index_ct[end+1][j][k]] = 0.0;
			dudx[par->index_ct[end+2][j][k]] = 0.0;
		}
	}
    
    
	return;
}


/* ******************************************************************* */
/* take a derivative in x (stream direction) for v (index_ved)		   */
/* derivatives are defined at corner (index_cn)						   */
/* index_ved[]: [start-1->end][0->ny][0->nz-1]                         */
/* index_cn[]: [start->end][0->ny][0->nz-1]                            */
/* ******************************************************************* */
void derivative_x1_ved2cn__2(
							 double *v,
							 double *dvdx,
							 parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
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
/* take a derivative in x (stream direction) for v (index_ved)		   */
/* derivatives are defined at corner (index_cn)						   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]					   */
/* ******************************************************************* */
void derivative_x1_ved2cn__4(
                             double *v,
                             double *dvdx,
                             parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dx = par->dx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdx[par->index_cn[i][j][k]] \
                = (v[par->index_ved[i][j][k]] - v[par->index_ved[i-1][j][k]])/dx;
			}
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/* take a derivative in x (stream direction) for v (index_ved)		   */
/* derivatives are defined at corner (index_cn)						   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]					   */
/* ******************************************************************* */
void derivative_x3_ved2cn__4(
                             double *v,
                             double *dvdx,
                             parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dx = par->dx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdx[par->index_cn[i][j][k]] \
                = (v[par->index_ved[i+1][j][k]] - v[par->index_ved[i-2][j][k]])/(3.0*dx);
			}
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/* take a derivative in x (stream direction) for w (index_ct)		   */
/* derivatives are defined at index_ued								   */
/* dwdx is not defined at i = start-1								   */
/* index_ct[]: [start-1->end][-1->ny][0->nz-1]					   */
/* index_ued[]: [start-1->end][-1->ny][0->nz-1]					   */
/* dwdx is not defined at i = start-1								   */
/* ******************************************************************* */
void derivative_x1_ct2ued__2(
							 double *w,
							 double *dwdx1,
							 parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dx = par->dx;
	
	for (i=start; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdx1[par->index_ued[i][j][k]] \
				= (w[par->index_ct[i][j][k]] - w[par->index_ct[i-1][j][k]])/dx;
			}
		}
	}
	
	i = start-1;
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			dwdx1[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	
	return;
}

/* ******************************************************************* */
/* take a derivative in x (stream direction) for w (index_ct)		   */
/* derivatives are defined at index_ued								   */
/* dwdx is not defined at i = start-1								   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* dwdx is not defined at i = start-3								   */
/* ******************************************************************* */
void derivative_x1_ct2ued__4(
                             double *w,
                             double *dwdx1,
                             parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dx = par->dx;
	
	for (i=start-2; i<=end+2; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdx1[par->index_ued[i][j][k]] \
                = (w[par->index_ct[i][j][k]] - w[par->index_ct[i-1][j][k]])/dx;
			}
		}
	}
	
	i = start-3;
	for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			dwdx1[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
    
	return;
}

/* ******************************************************************* */
/* take a derivative in x (stream direction) for w (index_ct)		   */
/* derivatives are defined at index_ued								   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* dwdx is not defined at i = start-3, start-2, end+2				   */
/* ******************************************************************* */
void derivative_x3_ct2ued__4(
                             double *w,
                             double *dwdx3,
                             parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dx = par->dx;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdx3[par->index_ued[i][j][k]] \
                = (w[par->index_ct[i+1][j][k]] - w[par->index_ct[i-2][j][k]])/(3.0*dx);
			}
		}
	}
	
    for (j=-2; j<=ny+1; ++j){
		for (k=0; k<=nz-1; ++k){
			dwdx3[par->index_ued[start-3][j][k]] = 0.0;
			dwdx3[par->index_ued[start-2][j][k]] = 0.0;
			dwdx3[par->index_ued[end+2][j][k]] = 0.0;
		}
	}
	
    
	return;
}


/* ******************************************************************* */
/* take a derivative in y (normal direction) for u (index_ued)		   */
/* derivatives are defined at corner (index_cn)						   */
/* index_ued[]: [start-1->end][-1->ny][0->nz-1]					   */
/* index_cn[]: [start->end][0->ny][0->nz-1]					   */
/* ******************************************************************* */
void derivative_y1_ued2cn__2(
							 double *u,
							 double *dudy,
							 parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
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
/* take a derivative in y (normal direction) for u (index_ued)		   */
/* derivatives are defined at corner (index_cn)						   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]					   */
/* ******************************************************************* */
void derivative_y1_ued2cn__4(
                             double *u,
                             double *dudy,
                             parameters *par)
{
    
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dudy[par->index_cn[i][j][k]] \
                = (u[par->index_ued[i][j][k]] - u[par->index_ued[i][j-1][k]])/dy;
			}
		}
	}
    
	
	return;
}


/* ******************************************************************* */
/* take a derivative in y (normal direction) for u (index_ued)		   */
/* derivatives are defined at corner (index_cn)						   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_cn[]: [start-1->end+1][-1->ny+1][0->nz-1]					   */
/* u_cn is not defined at j = -1, ny+1								   */
/* ******************************************************************* */
void derivative_y3_ued2cn__4(
                             double *u,
                             double *dudy,
                             parameters *par)
{
    
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-1; i<=end+1; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dudy[par->index_cn[i][j][k]] \
                = (u[par->index_ued[i][j+1][k]] - u[par->index_ued[i][j-2][k]])/(3.0*dy);
			}
		}
	}
	
	//j = -1;
	//j = ny+1;
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			dudy[par->index_cn[i][-1][k]] = 0.0;
			dudy[par->index_cn[i][ny+1][k]] = 0.0;
		}
	}
	
	return;
}


/* ******************************************************************* */
/* take a derivative in y (normal direction) for v (index_ved)		   */
/* derivatives are defined at center, index_ct  					   */
/* index_ved[]: [start-1->end][0->ny][0->nz-1]					   */
/* index_ct[]: [start-1->end][-1->ny][0->nz-1]					   */
/* at j = -1, ny, values are not defined								*/
/* ******************************************************************* */
void derivative_y1_ved2ct__2(
							 double *v,
							 double *dvdy,
							 parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdy[par->index_ct[i][j][k]] \
				= (v[par->index_ved[i][j+1][k]] - v[par->index_ved[i][j][k]])/dy;
			}
		}
	}
	
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -1;
			dvdy[par->index_ct[i][j][k]] = 0.0;
			j = ny;
			dvdy[par->index_ct[i][j][k]] = 0.0;
		}
	}
    
	
	
	return;
}


/* ******************************************************************* */
/* take a derivative in y (normal direction) for v (index_ved)		   */
/* derivatives are defined at center, index_ct  					   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* ******************************************************************* */
void derivative_y1_ved2ct__4(
                             double *v,
                             double *dvdy,
                             parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-3; i<=end+2; ++i){
		for (j=-2; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdy[par->index_ct[i][j][k]] \
                = (v[par->index_ved[i][j+1][k]] - v[par->index_ved[i][j][k]])/dy;
			}
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/* take a derivative in y (normal direction) for v (index_ved)		   */
/* derivatives are defined at center, index_ct  					   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* v_ct is not defined at j = -2, j = ny+1							   */
/* ******************************************************************* */
void derivative_y3_ved2ct__4(
                             double *v,
                             double *dvdy,
                             parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-3; i<=end+2; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dvdy[par->index_ct[i][j][k]] \
                = (v[par->index_ved[i][j+2][k]] - v[par->index_ved[i][j-1][k]])/(3.0*dy);
			}
		}
	}
	
    for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			dvdy[par->index_ct[i][-2][k]] = 0.0;
			dvdy[par->index_ct[i][ny+1][k]] = 0.0;
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/* take a derivative in y (normal direction) for w (index_ct)		   */
/* derivatives are defined at v edge, index_ved						   */
/* index_ct[]: [start-1->end][-1->ny][0->nz-1]					   */
/* index_ved[]: [start-1->end][0->ny][0->nz-1]						   */
/* ******************************************************************* */
void derivative_y1_ct2ved__2(
							 double *w,
							 double *dwdy1,
							 parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdy1[par->index_ved[i][j][k]] \
				= (w[par->index_ct[i][j][k]] - w[par->index_ct[i][j-1][k]])/dy;
			}
		}
	}
	
	return;
}




/* ******************************************************************* */
/* take a derivative in y (normal direction) for w (index_ct)		   */
/* derivatives are defined at v edge, index_ved						   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* w_ved is not defined at j = -2, ny+2								   */
/* ******************************************************************* */
void derivative_y1_ct2ved__4(
                             double *w,
                             double *dwdy1,
                             parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-3; i<=end+2; ++i){
		for (j=-1; j<=ny+1; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdy1[par->index_ved[i][j][k]] \
                = (w[par->index_ct[i][j][k]] - w[par->index_ct[i][j-1][k]])/dy;
			}
		}
	}
	
    for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			dwdy1[par->index_ved[i][-2][k]] = 0.0;
			dwdy1[par->index_ved[i][ny+2][k]] = 0.0;
		}
	}
	
	return;
}


/* ******************************************************************* */
/* take a derivative in y (normal direction) for w (index_ct)		   */
/* derivatives are defined at v edge, index_ved						   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* w_ved is not defined at j = -2, -1, ny+1, ny+2					   */
/* ******************************************************************* */
void derivative_y3_ct2ved__4(
                             double *w,
                             double *dwdy3,
                             parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	const double dy = par->dy;
	
	for (i=start-3; i<=end+2; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				dwdy3[par->index_ved[i][j][k]] \
                = (w[par->index_ct[i][j+1][k]] - w[par->index_ct[i][j-2][k]])/(3.0*dy);
			}
		}
	}
	
	
	for (i=start-3; i<=end+2; ++i){
		for (k=0; k<=nz-1; ++k){
			dwdy3[par->index_ved[i][-2][k]] = 0.0;
			dwdy3[par->index_ved[i][-1][k]] = 0.0;
            dwdy3[par->index_ved[i][ny+1][k]] = 0.0;
			dwdy3[par->index_ved[i][ny+2][k]] = 0.0;
		}
	}
	
    
	return;
}




/* ******************************************************************* */
/* take a derivative in z (fft direction) for u (defined at ued)	   */
/* derivatives are also defined at u edge, index_ued				   */
/* index_ued[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* ******************************************************************* */
void derivative_z_ued2ued(
                          double *u,
                          double *dudz,
                          parameters *par)
{
    
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
    
	
	for (i=start-3; i<=end+2; ++i){
		for (j=-2; j<=ny+1; ++j){
			k = 0;
			dudz[par->index_ued[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				dudz[par->index_ued[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				dudz[par->index_ued[i][j][nz-k]] = u[par->index_ued[i][j][k]]*par->kz[k];
				dudz[par->index_ued[i][j][k]] =  - u[par->index_ued[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/* take a derivative in z (fft direction) for u (defined at ued)	   */
/* derivatives are also defined at u edge, index_ued				   */
/* index_ued[]: [start-1->end][-1->ny][0->nz-1]					   */
/* ******************************************************************* */
void derivative_z_ued2ued__2(
                             double *u,
                             double *dudz,
                             parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	
	
	for (i=start-1; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			k = 0;
			dudz[par->index_ued[i][j][k]] = 0.0;
			
			if (nz%2 == 0){
				k = nz/2;
				dudz[par->index_ued[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				dudz[par->index_ued[i][j][nz-k]] = u[par->index_ued[i][j][k]]*par->kz[k];
				dudz[par->index_ued[i][j][k]] =  - u[par->index_ued[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/* take a derivative in z (fft direction) for v (defined at ved)	   */
/* derivatives are also defined at v edge, index_ved				   */
/* index_ved[]: [start-3->end+2][-2->ny+2][0->nz-1]					   */
/* ******************************************************************* */
void derivative_z_ved2ved(
                          double *v,
                          double *dvdz,
                          parameters *par)
{
    
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
    
	for (i=start-3; i<=end+2; ++i){
		for (j=-2; j<=ny+2; ++j){
			k = 0;
			dvdz[par->index_ved[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				dvdz[par->index_ved[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				dvdz[par->index_ved[i][j][nz-k]] = v[par->index_ved[i][j][k]]*par->kz[k];
				dvdz[par->index_ved[i][j][k]] =  - v[par->index_ved[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/* take a derivative in z (fft direction) for v (defined at ved)	   */
/* derivatives are also defined at v edge, index_ved				   */
/* index_ved[]: [start-1->end][0->ny][0->nz-1]					   */
/* ******************************************************************* */
void derivative_z_ved2ved__2(
                             double *v,
                             double *dvdz,
                             parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			k = 0;
			dvdz[par->index_ved[i][j][k]] = 0.0;
			
			if (nz%2 == 0){
				k = nz/2;
				dvdz[par->index_ved[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				dvdz[par->index_ved[i][j][nz-k]] = v[par->index_ved[i][j][k]]*par->kz[k];
				dvdz[par->index_ved[i][j][k]] =  - v[par->index_ved[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
	return;
}

/* ******************************************************************* */
/* take a derivative in z (fft direction) for w (defined at ct)	       */
/* derivatives are also defined at center, index_ct					   */
/* index_ct[]: [start-3->end+2][-2->ny+1][0->nz-1]					   */
/* ******************************************************************* */
void derivative_z_ct2ct(
                        double *w,
                        double *dwdz,
                        parameters *par)
{
    
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	
	for (i=start-3; i<=end+2; ++i){
		for (j=-2; j<=ny+1; ++j){
			k = 0;
			dwdz[par->index_ct[i][j][k]] = 0.0;
            
			if (nz%2 == 0){
				k = nz/2;
				dwdz[par->index_ct[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				dwdz[par->index_ct[i][j][nz-k]] = w[par->index_ct[i][j][k]]*par->kz[k];
				dwdz[par->index_ct[i][j][k]] =  - w[par->index_ct[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
	return;
}



/* ******************************************************************* */
/* take a derivative in z (fft direction) for w (defined at ct)	       */
/* derivatives are also defined at center, index_ct					   */
/* index_ct[]: [start-1->end][-1->ny][0->nz-1]					   */
/* ******************************************************************* */
void derivative_z_ct2ct__2(
                           double *w,
                           double *dwdz,
                           parameters *par)
{
	
	int i, j, k;
    const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	
	for (i=start-1; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			k = 0;
			dwdz[par->index_ct[i][j][k]] = 0.0;
			
			if (nz%2 == 0){
				k = nz/2;
				dwdz[par->index_ct[i][j][k]] = 0.0;
			}
			
			for (k=1; k<=(nz-1)/2; ++k){
				dwdz[par->index_ct[i][j][nz-k]] = w[par->index_ct[i][j][k]]*par->kz[k];
				dwdz[par->index_ct[i][j][k]] =  - w[par->index_ct[i][j][nz-k]]*par->kz[nz-k];
			}
		}
	}
	
	return;
}
