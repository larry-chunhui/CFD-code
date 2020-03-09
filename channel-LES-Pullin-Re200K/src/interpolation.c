/*
 *  interpolation.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Allocate and deallocate the memory for inter-> variables.
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

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


void init_interpolation(
                        interpolated *inter,
                        parameters *par)
{
	
	const int local_nx = par->local_nx;
	const int nz = par->nz;
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	const int local_size_les = par->local_size_les;
	
	
	inter->u_ued = FFTW_MALLOC(local_size_ued);  // freed
	inter->u_phys_ued = FFTW_MALLOC(local_size_ued);  // freed
    
	inter->v_ved = FFTW_MALLOC(local_size_ved);  // freed
	inter->v_phys_ved = FFTW_MALLOC(local_size_ved);  // freed
	
	inter->w_ct = FFTW_MALLOC(local_size_ct);  // freed
	inter->w_phys_ct = FFTW_MALLOC(local_size_ct);  // freed
    
	inter->u_phys_cn1 = FFTW_MALLOC(local_size_cn); // freed
	inter->v_phys_cn1 = FFTW_MALLOC(local_size_cn); // freed
	
	inter->u_phys_ct1 = FFTW_MALLOC(local_size_ct); // freed
	inter->v_phys_ct1 = FFTW_MALLOC(local_size_ct); // freed
	
	inter->w_phys_ved1 = FFTW_MALLOC(local_size_ved); // freed
	inter->w_phys_ued1 = FFTW_MALLOC(local_size_ued); // freed
	
	
	if (par->fd_order == 4){
		inter->u_phys_cn3 = FFTW_MALLOC(local_size_cn); // freed
		inter->v_phys_cn3 = FFTW_MALLOC(local_size_cn); // freed
		
		inter->u_phys_ct3 = FFTW_MALLOC(local_size_ct); // freed
		inter->v_phys_ct3 = FFTW_MALLOC(local_size_ct); // freed
		
		inter->w_phys_ved3 = FFTW_MALLOC(local_size_ved); // freed
		inter->w_phys_ued3 = FFTW_MALLOC(local_size_ued); // freed
		
		inter->uv3_top = FFTW_MALLOC((local_nx+3)*nz); // freed
		inter->uv3_bot = FFTW_MALLOC((local_nx+3)*nz); // freed
		inter->wv3_top = FFTW_MALLOC((local_nx+3)*nz); // freed
		inter->wv3_bot = FFTW_MALLOC((local_nx+3)*nz); // freed
	}
	
	
	
	inter->Txx_ct = FFTW_MALLOC(local_size_ct); // freed
	inter->Tyy_ct = FFTW_MALLOC(local_size_ct); // freed
	inter->Tzz_ct = FFTW_MALLOC(local_size_ct); // freed
	
	inter->Txy_cn = FFTW_MALLOC(local_size_cn); // freed
	inter->Tzx_ued = FFTW_MALLOC(local_size_ued); // freed
	inter->Tyz_ved = FFTW_MALLOC(local_size_ved); // freed
	
	
	if (par->les == 1){
		inter->u_phys_les = FFTW_MALLOC(local_size_les);//freed
		inter->v_phys_les = FFTW_MALLOC(local_size_les);//freed
		inter->w_phys_les = FFTW_MALLOC(local_size_les);//freed
	}
	
	if (par->les == 1 && par->bc == 1){
		inter->bc_duudx_t = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_duwdz_t = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_duudx_b = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_duwdz_b = FFTW_MALLOC((local_nx+1)*(nz));//freed
        
		inter->bc_ududx_t = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_wdudz_t = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_ududx_b = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_wdudz_b = FFTW_MALLOC((local_nx+1)*(nz));//freed
		
		inter->bc_uv_b = FFTW_MALLOC((local_nx+1)*(nz));//freed
		inter->bc_uv_t = FFTW_MALLOC((local_nx+1)*(nz));//freed
	}
	
    
    
	
	
	return;
}


void finalize_interpolation(
                            interpolated *inter,
                            parameters *par)
{
	
	fftw_free(inter->u_ued);
	fftw_free(inter->u_phys_ued);
	
	fftw_free(inter->v_ved);
	fftw_free(inter->v_phys_ved);
	
	fftw_free(inter->w_ct);
	fftw_free(inter->w_phys_ct);
	
	fftw_free(inter->u_phys_cn1);
	fftw_free(inter->u_phys_ct1);
    
	fftw_free(inter->v_phys_cn1);
	fftw_free(inter->v_phys_ct1);
	
	fftw_free(inter->w_phys_ued1);
	fftw_free(inter->w_phys_ved1);
	
	fftw_free(inter->Txx_ct);
	fftw_free(inter->Tyy_ct);
	fftw_free(inter->Tzz_ct);
	
	fftw_free(inter->Txy_cn);
	fftw_free(inter->Tzx_ued);
	fftw_free(inter->Tyz_ved);
	
	
	if (par->les == 1){
		fftw_free(inter->u_phys_les);
		fftw_free(inter->v_phys_les);
		fftw_free(inter->w_phys_les);
	}
    
	if (par->les == 1 && par->bc == 1){
		fftw_free(inter->bc_duudx_t);
		fftw_free(inter->bc_duwdz_t);
		fftw_free(inter->bc_duudx_b);
		fftw_free(inter->bc_duwdz_b);
        
		fftw_free(inter->bc_ududx_t);
		fftw_free(inter->bc_wdudz_t);
		fftw_free(inter->bc_ududx_b);
		fftw_free(inter->bc_wdudz_b);
		
		fftw_free(inter->bc_uv_t);
		fftw_free(inter->bc_uv_b);
	}
	
	
	
	
	if (par->fd_order == 4){
		fftw_free(inter->u_phys_cn3);
		fftw_free(inter->u_phys_ct3);
        
		fftw_free(inter->v_phys_cn3);
		fftw_free(inter->v_phys_ct3);
        
		fftw_free(inter->w_phys_ued3);
		fftw_free(inter->w_phys_ved3);
		
		fftw_free(inter->uv3_top);
		fftw_free(inter->uv3_bot);
		
		fftw_free(inter->wv3_top);
		fftw_free(inter->wv3_bot);
	}
	
	return;
}

