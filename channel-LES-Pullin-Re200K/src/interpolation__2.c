/*
 *  interpolation__2.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Second interpolation functions.
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

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


void interpolation__2(
                      int flag,
                      variables *var,
                      neighbors *shared,
                      interpolated *inter,
                      parameters *par,
                      fftwplans *planptr)
{
	
	int i, j, k;
	int index1, index2;
	double *copy_ct;
	double *copy_ued;
	double *copy_ved;
	
	bigger_array_u_ued__2(flag, var, shared, inter->u_ued, par);
	bigger_array_v_ved__2(flag, var, shared, inter->v_ved, par);
	bigger_array_wp_ct__2(flag, var, shared, inter->w_ct, par);
	
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
	
	
	interpolate1_ued_cn__2(inter->u_phys_ued, inter->u_phys_cn1, par);
	interpolate1_ued_ct__2(inter->u_phys_ued, inter->u_phys_ct1, par);
	
	interpolate1_ved_cn__2(inter->v_phys_ved, inter->v_phys_cn1, par);
	interpolate1_ved_ct__2(inter->v_phys_ved, inter->v_phys_ct1, par);
	
	interpolate1_ct_ued__2(inter->w_phys_ct, inter->w_phys_ued1, par);
	interpolate1_ct_ved__2(inter->w_phys_ct, inter->w_phys_ved1, par);
	
	
	
	double *w_phys_cn;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	if (par->les == 1){
		
		w_phys_cn = FFTW_MALLOC(par->local_size_cn); //freed
		
		// i = start is not defined
		interpolate_ct2cn__2(inter->w_phys_ct, w_phys_cn, par);
		
		for (i=start; i<=end; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					index1 = par->index_les[i][j][k];
					index2 = par->index_cn[i][j][k];
					
					inter->u_phys_les[index1] = inter->u_phys_cn1[index2];
					
					inter->v_phys_les[index1] = inter->v_phys_cn1[index2];
					
					inter->w_phys_les[index1] = w_phys_cn[index2];
					
				}
			}
		}
		
		// values at i = start-1 need to be obtained from neighboring processes
		{
			int index;
			double *in_l;
			
			in_l = dvector(0, 1*(ny+1)*nz - 1);
			
			share_les_periodic__2(par->my_rank, inter->u_phys_les, in_l, par, 401);
			i = start-1;
			index = 0;
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					inter->u_phys_les[par->index_les[i][j][k]] = in_l[index];
					index += 1;
				}
			}
			
			
			share_les_periodic__2(par->my_rank, inter->v_phys_les, in_l, par, 501);
			
			i = start-1;
			index = 0;
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					inter->v_phys_les[par->index_les[i][j][k]] = in_l[index];
					index += 1;
				}
			}
			
			
			share_les_periodic__2(par->my_rank, inter->w_phys_les, in_l, par, 601);
			
			i = start-1;
			index = 0;
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					inter->w_phys_les[par->index_les[i][j][k]] = in_l[index];
					index += 1;
				}
			}
			
			
			free_dvector(in_l, 0, 1*(ny+1)*nz-1);
		}
		
		
		fftw_free(w_phys_cn);
	}
    
	return;
	
	
}


/* ******************************************************************* */
/*
 rearrange velocity w[]
 index_ct[]: [start-1->end][-1->ny][0->nz-1]
 index_ved[]: [start-1->end][0->ny][0->nz-1]
 
 w_ved is not defined
 at j = 0, ny
 */
/* ******************************************************************* */
void interpolate1_ct_ved__2(
                            double *w_ct,
                            double *w_ved,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				w_ved[par->index_ved[i][j][k]] \
                = (w_ct[par->index_ct[i][j-1][k]] + w_ct[par->index_ct[i][j][k]])/2.0;
			}
		}
	}
	
    
	j = 0;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			w_ved[par->index_ved[i][j][k]] = 0.0;
		}
	}
	
    
	j = ny;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			w_ved[par->index_ved[i][j][k]] = 0.0;
		}
	}
	
	return;
}




/* ******************************************************************* */
/*
 rearrange velocity w[]
 index_ct[]: [start-1->end][-1->ny][0->nz-1]
 index_ued[]: [start-1->end][-1->ny][0->nz-1]
 
 w_ued is not defined
 at i = start-1
 */
/* ******************************************************************* */
void interpolate1_ct_ued__2(
                            double *w_ct,
                            double *w_ued,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				w_ued[par->index_ued[i][j][k]] \
                = (w_ct[par->index_ct[i-1][j][k]] + w_ct[par->index_ct[i][j][k]])/2.0;
			}
		}
	}
	
	
	i = start-1;
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			w_ued[par->index_ued[i][j][k]] = 0.0;
		}
	}
    
	return;
}



/* ******************************************************************* */
/*
 rearrange velocity v[]
 index_ved[]: [start-1->end][0->ny][0->nz-1]
 index_cn[]: [start->end][0->ny][0->nz-1]
 */
/* ******************************************************************* */
void interpolate1_ved_cn__2(
                            double *v_ved,
                            double *v_cn,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end; ++i){
		for (j=0; j<=ny; ++j){
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
 rearrange velocity v[] (index_v[][][])
 index_ved[]: [start-1->end][0->ny][0->nz-1]
 index_ct[]: [start-1->end][-1->ny][0->nz-1]
 
 v_ct is not defined at j = -1, ny
 */
/* ******************************************************************* */
void interpolate1_ved_ct__2(
                            double *v_ved,
                            double *v_ct,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				v_ct[par->index_ct[i][j][k]] \
                = (v_ved[par->index_ved[i][j][k]] + v_ved[par->index_ved[i][j+1][k]])/2.0;
			}
		}
	}
	
	j = -1;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			v_ct[par->index_ct[i][j][k]] = 0.0;
		}
	}
	
	j = ny;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			v_ct[par->index_ct[i][j][k]] = 0.0;
		}
	}
	
    
	
	
    
	return;
}



/* ******************************************************************* */
/*
 rearrange velocity u[]
 index_ued[]: [-1->nx][-1->ny][0->nz-1]
 index_cn[]: [start->end][0->ny][0->nz-1]
 */
/* ******************************************************************* */
void interpolate1_ued_cn__2(
                            double *u_ued,
                            double *u_cn,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end; ++i){
		for (j=0; j<=ny; ++j){
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
 rearrange velocity u[] (index_u[][][])
 index_ued[]: [start-1->end][-1->ny][0->nz-1]
 index_ct[]: [start-1->end][-1->ny][0->nz-1]
 u_ct is not defined at i = end
 */
/* ******************************************************************* */
void interpolate1_ued_ct__2(
                            double *u_ued,
                            double *u_ct,
                            parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    
	for (i=start-1; i<=end-1; ++i){
		for (j=-1; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				u_ct[par->index_ct[i][j][k]] \
                = (u_ued[par->index_ued[i][j][k]] + u_ued[par->index_ued[i+1][j][k]])/2.0;
			}
		}
	}
	
	i = end;
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			u_ct[par->index_ct[i][j][k]] = 0.0;
		}
	}
    
	return;
}




/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start->end][0->ny][0->nz-1]
 index_ued[]: [start-1->end][-1->ny][0->nz-1]
 
 i = start-1 will not be defined
 j = -1, ny will not be defined
 */
/* ******************************************************************* */
void interpolate_y1_cn2ued__2(
							  double *u_cn,
							  double *u_ed,
							  parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	for (i=start; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				u_ed[par->index_ued[i][j][k]] \
				= (u_cn[par->index_cn[i][j][k]] + u_cn[par->index_cn[i][j+1][k]])/2.0;
			}
		}
	}
	
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			i = start-1;
			u_ed[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			j = -1;
			u_ed[par->index_ued[i][j][k]] = 0.0;
			j = ny;
			u_ed[par->index_ued[i][j][k]] = 0.0;
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/*
 interpolate velocity
 index_cn[]: [start->end][0->ny][0->nz-1]
 index_ved[]: [start-1->end][0->ny][0->nz-1]
 
 ved is not defined
 at i = start-1, end
 */
/* ******************************************************************* */
void interpolate_x1_cn2ved__2(
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
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				ved[par->index_ved[i][j][k]] \
				= (v_cn[par->index_cn[i][j][k]] + v_cn[par->index_cn[i+1][j][k]])/2.0;
			}
		}
	}
	
	for (j=0; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			i = start-1;
			ved[par->index_ved[i][j][k]] = 0.0;
			i = end;
			ved[par->index_ved[i][j][k]] = 0.0;
		}
	}
    
	
	return;
}




/* ******************************************************************* */
/*
 rearrange velocity u[] (index_u[][][])
 index_ued[]: [start-1->end][-1->ny][0->nz-1]
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_x1_ued2wp__2(
                              double *u_ued,
                              double *u_wp,
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
				u_wp[par->index_wp[i][j][k]] \
				= (u_ued[par->index_ued[i][j][k]] + u_ued[par->index_ued[i+1][j][k]])/2.0;
			}
		}
	}
	
	return;
}



/* ******************************************************************* */
/*
 by interpolate1_ct_ved__2()
 index_ct[]: [start-1->end][-1->ny][0->nz-1]
 index_ved[]: [start-1->end][0->ny][0->nz-1]
 
 ved is not defined
 at j = 0, ny
 
 interpolate1_ved_cn__2()
 index_ved[]: [start-1->end][0->ny][0->nz-1]
 index_cn[]: [start->end][0->ny][0->nz-1]
 
 cn is not defined
 at j = 0, ny
 */


/*
 index_ct[]: [start-1->end][-1->ny][0->nz-1]
 index_ued[]: [start-1->end][-1->ny][0->nz-1]
 
 w_ued is not defined
 at i = start-1
 
 index_ued[]: [start-1->end][-1->ny][0->nz-1]
 index_cn[]: [start->end][0->ny][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_ct2cn__2(
						  double *ct,
						  double *cn,
						  parameters *par)
{
	const int local_size_ued = par->local_size_ued;
	
	double *scratch;
	
	// interpolate ct data onto u_edge
	scratch = FFTW_MALLOC(local_size_ued);
	interpolate1_ct_ued__2(ct, scratch, par);
	
	// further interpolate onto corner
	interpolate1_ued_cn__2(scratch, cn, par);
	
	fftw_free(scratch);
    
    
	
	return;
}



/* ******************************************************************* */
/*
 index_wp[]: [start->end-1][0->ny-1][0->nz-1]
 */
/* ******************************************************************* */
void interpolate_cn2wp__2(
						  double *cn,
						  double *wp,
						  parameters *par)
{
	const int local_size_ued = par->local_size_ued;
	
	double *scratch1;
	
	// interpolate dudy onto u_edge ******************
	//index_cn[]: [start->end][0->ny][0->nz-1]
	//index_ued[]: [start-1->end][-1->ny][0->nz-1]
	//i = start-1 will not be defined
	//j = -1, ny will not be defined
	
	scratch1 = FFTW_MALLOC(local_size_ued);
	
	interpolate_y1_cn2ued__2(cn, scratch1, par);
	
	
	
	// further interpolate onto center ******************
	// For scratch1,
	// index_ued[]: [start-1->end][-1->ny][0->nz-1]
	// i = start will not be defined
	// to define interpolated wp
	
	interpolate_x1_ued2wp__2(scratch1, wp, par);
	
	fftw_free(scratch1);
	
	
	return;
}


/* ******************************************************************* */
/*
 rearrange velocity u[] (index_u[][][])
 
 index_ued[] is defined at the right edge of each cell similar to u[], but boundary
 values are included.
 index_ued[]: [-1->nx][-1->ny][0->nz-1]
 
 flag = 0 (homogeneous boundary)
 flag = 1 (inhomogeneous boudary)
 */
/* ******************************************************************* */
void bigger_array_u_ued__2(
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
	
    
	double *u_r;
	double *u_l;
    
	u_r = dvector(0, (ny+2)*nz-1);
	u_l = dvector(0, (ny+2)*nz-1);
	
	for (i=0; i<=(ny+2)*nz-1; ++i){
		u_r[i] = 0.0;
		u_l[i] = 0.0;
	}
	
	
	index = 0;
	// bottom ghost cell
	for (k=0; k<=nz-1; ++k){
		u_r[index] = 2.0*flag*shared->u_in_r[index] \
        - shared->u_in_r[index+nz];
		index = index + 1;
	}
	// interior
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			u_r[index] = shared->u_in_r[index];
			index = index + 1;
		}
	}
	// top ghost cell
	for (k=0; k<=nz-1; ++k){
		u_r[index] = 2.0*flag*shared->u_in_r[index] \
        - shared->u_in_r[index-nz];
		index = index + 1;
	}
	
	
	index = 0;
	// bottom ghost cell
	for (k=0; k<=nz-1; ++k){
		u_l[index] = 2.0*flag*shared->u_in_l[index] \
        - shared->u_in_l[index+nz];
		index = index + 1;
	}
	// interior
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			u_l[index] = shared->u_in_l[index];
			index = index + 1;
		}
	}
	// top ghost cell
	for (k=0; k<=nz-1; ++k){
		u_l[index] = 2.0*flag*shared->u_in_l[index] \
        - shared->u_in_l[index-nz];
		index = index + 1;
	}
	
    
    
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// index_ued[]: [-1->nx][-1->ny][0->nz-1]
	// define ued at the edges
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				// left boundary
				if (i == start-1){
					u_ed[par->index_ued[i][j][k]] = u_l[nz*(j+1) + k];
                    // right boundary
				}else if (i == end){
					u_ed[par->index_ued[i][j][k]] = u_r[nz*(j+1) + k];
                    // interior
				}else{
					u_ed[par->index_ued[i][j][k]] = var->u[par->index_u[i][j][k]];
				}
			}
		}
	}
	
	// bottom boundary
	j = -1;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i == start-1){
				u_ed[par->index_ued[i][j][k]] = u_l[nz*(j+1) + k];
                // right boundary
			}else if (i == end){
				u_ed[par->index_ued[i][j][k]] = u_r[nz*(j+1) + k];
                // interior
			}else{
				u_ed[par->index_ued[i][j][k]] = 2.0*flag*var->bc_u_b[par->index_tb[i][k]] \
                - var->u[par->index_u[i][0][k]];
			}
		}
	}
	
	
	
	// top boundary
	j = ny;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i == start-1){
				u_ed[par->index_ued[i][j][k]] = u_l[nz*(j+1) + k];
                // right boundary
			}else if (i == end){
				u_ed[par->index_ued[i][j][k]] = u_r[nz*(j+1) + k];
                // interior
			}else{
				u_ed[par->index_ued[i][j][k]] = 2.0*flag*var->bc_u_t[par->index_tb[i][k]] \
                - var->u[par->index_u[i][ny-1][k]];
			}
		}
	}
	
	
	free_dvector(u_r, 0, (ny+2)*nz-1);
	free_dvector(u_l, 0, (ny+2)*nz-1);
    
	return;
}


/* ******************************************************************* */
/*
 rearrange velocity v[] (index_v[][][])
 
 index_ved[] is defined at the top edge of each cell similar to v[], but boundary
 values are included. index_ved[]: [-1->nx][0->ny][0->nz-1]
 
 flag = 0 (homogeneous boundary)
 flag = 1 (inhomogeneous boudary)
 */
/* ******************************************************************* */
void bigger_array_v_ved__2(
                           double flag,
                           variables *var,
                           neighbors *shared,
                           double *v_ed,
                           parameters *par)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *v_r;
	double *v_l;
    
	v_r = dvector(0, (ny+1)*nz-1);
	v_l = dvector(0, (ny+1)*nz-1);
	
	for (i=0; i<=(ny+1)*nz-1; ++i){
		v_r[i] = 0.0;
		v_l[i] = 0.0;
	}
	
	
	index = 0;
	for (j=0; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			if (flag == 0.0 && (j == 0 || j == ny)){
				v_r[index] = 0.0;
			}else{
				v_r[index] = shared->v_in_r[index];
			}
			index = index + 1;
		}
	}
	
	index = 0;
	for (j=0; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			if (flag == 0.0 && (j == 0 || j == ny)){
				v_l[index] = 0.0;
			}else{
				v_l[index] = shared->v_in_l[index];
			}
			index = index + 1;
		}
	}
    
	
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// index_ved[]: [-1->nx][0->ny][0->nz-1]
	// define ved at the edges
	for (i=start-1; i<=end; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				// left boundary
				if (i == start-1){
					v_ed[par->index_ved[i][j][k]] = v_l[nz*j + k];
                    // right boundary
				}else if (i == end){
					v_ed[par->index_ved[i][j][k]] = v_r[nz*j + k];
                    // interior
				}else{
					v_ed[par->index_ved[i][j][k]] = var->v[par->index_v[i][j][k]];
				}
			}
		}
	}
	
	// bottom boundary
	j = 0;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i == start-1){
				v_ed[par->index_ved[i][j][k]] = v_l[nz*j + k];
                // right boundary
			}else if (i == end){
				v_ed[par->index_ved[i][j][k]] = v_r[nz*j + k];
                // interior
			}else{
				v_ed[par->index_ved[i][j][k]] = flag*var->bc_v_b[par->index_tb[i][k]];
			}
		}
	}
    
	
	
	// top boundary
	j = ny;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i == start-1){
				v_ed[par->index_ved[i][j][k]] = v_l[nz*j + k];
                // right boundary
			}else if (i == end){
				v_ed[par->index_ved[i][j][k]] = v_r[nz*j + k];
                // interior
			}else{
				v_ed[par->index_ved[i][j][k]] = flag*var->bc_v_t[par->index_tb[i][k]];
			}
		}
	}
	
	
	free_dvector(v_r, 0, (ny+1)*nz-1);
	free_dvector(v_l, 0, (ny+1)*nz-1);
    
	return;
}



/* ******************************************************************* */
/*
 rearrange velocity w[] (index_wp[][][])
 
 index_ct[] is defined as w[], but boundary values are included.
 index_ct[]: [-1->nx][-1->ny][0->nz-1]
 
 flag = 0 (homogeneous boundary)
 flag = 1 (inhomogeneous boudary)
 */
/* ******************************************************************* */
void bigger_array_wp_ct__2(
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
	
	double *w_r;
	double *w_l;
    
	w_r = dvector(0, (ny+2)*nz-1);
	w_l = dvector(0, (ny+2)*nz-1);
	
	for (i=0; i<=(ny+2)*nz-1; ++i){
		w_r[i] = 0.0;
		w_l[i] = 0.0;
	}
	
    
	index = 0;
	// bottom ghost cell
	for (k=0; k<=nz-1; ++k){
		w_r[index] = 2.0*flag*shared->w_in_r[index] \
        - shared->w_in_r[index+nz];
		index = index + 1;
	}
	// interior
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			w_r[index] = shared->w_in_r[index];
			index = index + 1;
		}
	}
	// top ghost cell
	for (k=0; k<=nz-1; ++k){
		w_r[index] = 2.0*flag*shared->w_in_r[index] \
        - shared->w_in_r[index-nz];
		index = index + 1;
	}
    
	
	index = 0;
	// bottom ghost cell
	for (k=0; k<=nz-1; ++k){
		w_l[index] = 2.0*flag*shared->w_in_l[index] \
        - shared->w_in_l[index+nz];
		index = index + 1;
	}
	// interior
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			w_l[index] = shared->w_in_l[index];
			index = index + 1;
		}
	}
	// top ghost cell
	for (k=0; k<=nz-1; ++k){
		w_l[index] = 2.0*flag*shared->w_in_l[index] \
        - shared->w_in_l[index-nz];
		index = index + 1;
	}
    
    
    
    
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// index_ct[]: [-1->nx][-1->ny][0->nz-1]
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				// left boundary
				if (i == start-1){
					w_ct[par->index_ct[i][j][k]] = w_l[nz*(j+1) + k];
                    // right boundary
				}else if (i == end){
					w_ct[par->index_ct[i][j][k]] = w_r[nz*(j+1) + k];
                    // interior
				}else{
					w_ct[par->index_ct[i][j][k]] = var->w[par->index_wp[i][j][k]];
				}
			}
		}
	}
	
	// bottom boundary
	j = -1;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i == start-1){
				w_ct[par->index_ct[i][j][k]] = w_l[nz*(j+1) + k];
                // right boundary
			}else if (i == end){
				w_ct[par->index_ct[i][j][k]] = w_r[nz*(j+1) + k];
                // interior
			}else{
				w_ct[par->index_ct[i][j][k]] = 2.0*flag*var->bc_w_b[par->index_tb[i][k]] \
                - var->w[par->index_wp[i][j+1][k]];
			}
		}
	}
	
	
	
	// top boundary
	j = ny;
	for (i=start-1; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			// left boundary
			if (i == start-1){
				w_ct[par->index_ct[i][j][k]] = w_l[nz*(j+1) + k];
                // right boundary
			}else if (i == end){
				w_ct[par->index_ct[i][j][k]] = w_r[nz*(j+1) + k];
                // interior
			}else{
				w_ct[par->index_ct[i][j][k]] = 2.0*flag*var->bc_w_t[par->index_tb[i][k]] \
                - var->w[par->index_wp[i][j-1][k]];
			}
		}
	}
    
	
	
	
	free_dvector(w_r, 0, (ny+2)*nz-1);
	free_dvector(w_l, 0, (ny+2)*nz-1);
    
	return;
}
