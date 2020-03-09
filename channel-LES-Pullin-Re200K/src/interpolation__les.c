/*
 *  interpolation__les.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Interpolations associated with LES.
 * void interpolation_les(variables *var, interpolated *inter, parameters *par)
 * void interpolation_les__2(variables *var, interpolated *inter, parameters *par)
 * void interpolate_Tij_ct__4(double *Tij_cn, double *Tij_ct, parameters *par)
 * void interpolate_Tij_ct__2(double *Tij_cn, double *Tij_ct, parameters *par)
 * _fill functions for Tij.
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
#include "interpolation__les.h"
#include "interpolation__4.h"
#include "interpolation__2.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

void interpolation_les(variables *var, interpolated *inter, parameters *par)
{
    
	int i;
	const int local_size_cn = par->local_size_cn;
	
	double *scratch;
	
	// from neighboring process
	Tij_fill_cn(var->Txx, par, 200);
	Tij_fill_cn(var->Tyy, par, 201);
	Tij_fill_cn(var->Tzz, par, 202);
	
	Tij_fill_cn(var->Txy, par, 203);
	Tij_fill_cn(var->Tyz, par, 204);
	Tij_fill_cn(var->Tzx, par, 205);
	
	// Tij defined at corner; Txy
	for (i=0; i<=local_size_cn-1; ++i){
		inter->Txy_cn[i] = var->Txy[i];
	}
	
	// Tij defined at ued; Tzx
	scratch = FFTW_MALLOC(par->local_size_ued); //freed
	
	interpolate_y1_cn2ued__4(var->Tzx, scratch, par);
	interpolate_y3_cn2ued__4(var->Tzx, inter->Tzx_ued, par);
	
	for (i=0; i<=par->local_size_ued-1; ++i){
		inter->Tzx_ued[i] = (9.0*scratch[i] - 1.0*inter->Tzx_ued[i])/8.0;
	}
	
	fftw_free(scratch);
	
	// Tij defined at ved; Tyz
	scratch = FFTW_MALLOC(par->local_size_ved); //freed
	
	interpolate_x1_cn2ved__4(var->Tyz, scratch, par);
	interpolate_x3_cn2ved__4(var->Tyz, inter->Tyz_ved, par);
	
	for (i=0; i<=par->local_size_ved-1; ++i){
		inter->Tyz_ved[i] = (9.0*scratch[i] - 1.0*inter->Tyz_ved[i])/8.0;
	}
	
	fftw_free(scratch);
	
	// Tij defined at center; Txx, Tyy, Tzz
	interpolate_Tij_ct__4(var->Txx, inter->Txx_ct, par);
	interpolate_Tij_ct__4(var->Tyy, inter->Tyy_ct, par);
	interpolate_Tij_ct__4(var->Tzz, inter->Tzz_ct, par);
	
	Txx_fill_ct(inter->Txx_ct, par, 301);
    
	Tyy_fill_ct(var->Tyy, inter->Tyy_ct, par);
    
	return;
}


void interpolation_les__2(variables *var, interpolated *inter, parameters *par)
{
	
	int i;
	const int local_size_cn = par->local_size_cn;
	
	// from neighboring process
	Tij_fill_cn__2(var->Txx, par, 200);
	Tij_fill_cn__2(var->Tyy, par, 201);
	Tij_fill_cn__2(var->Tzz, par, 202);
	
	Tij_fill_cn__2(var->Txy, par, 203);
	Tij_fill_cn__2(var->Tyz, par, 204);
	Tij_fill_cn__2(var->Tzx, par, 205);
	
	// Tij defined at corner; Txy
	for (i=0; i<=local_size_cn-1; ++i){
		inter->Txy_cn[i] = var->Txy[i];
	}
	
	// Tij defined at ued; Tzx
	interpolate_y1_cn2ued__2(var->Tzx, inter->Tzx_ued, par);
	
	// Tij defined at ved; Tyz
	interpolate_x1_cn2ved__2(var->Tyz, inter->Tyz_ved, par);
	// so far only ved:[start->end-1][0->ny][0->nz-1] are filled
	
	
	// Tij defined at center; Txx, Tyy, Tzz
	interpolate_Tij_ct__2(var->Txx, inter->Txx_ct, par);
	interpolate_Tij_ct__2(var->Tyy, inter->Tyy_ct, par);
	interpolate_Tij_ct__2(var->Tzz, inter->Tzz_ct, par);
	
	
	Txx_fill_ct__2(inter->Txx_ct, par, 301);
	
	return;
}




void interpolate_Tij_ct__4(double *Tij_cn, double *Tij_ct, parameters *par)
{
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *scratch;
	
	scratch = FFTW_MALLOC(par->local_size_wp); //freed
    
	interpolate_cn2wp__4(Tij_cn, scratch, par);
	
	for (i=0; i<=par->local_size_ct-1; ++i){
		Tij_ct[i] = 0.0;
	}
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Tij_ct[par->index_ct[i][j][k]]
                = scratch[par->index_wp[i][j][k]];
			}
		}
	}
	
	
	
	fftw_free(scratch);
	
	
	return;
}



void interpolate_Tij_ct__2(double *Tij_cn, double *Tij_ct, parameters *par)
{
	
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *scratch;
	
	scratch = FFTW_MALLOC(par->local_size_wp); //freed
	
	interpolate_cn2wp__2(Tij_cn, scratch, par);
	
	for (i=0; i<=par->local_size_ct-1; ++i){
		Tij_ct[i] = 0.0;
	}
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Tij_ct[par->index_ct[i][j][k]]
				= scratch[par->index_wp[i][j][k]];
			}
		}
	}
	
	fftw_free(scratch);
	
	
	return;
}

void Tij_fill_cn(double *Tij, parameters *par, int tag)
{
	
	int i, j, k;
	int index;
	int ijk;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *in_r, *in_l;
	
	in_l = dvector(0, 1*(ny+1)*nz - 1); //freed
	in_r = dvector(0, 2*(ny+1)*nz - 1); //freed
	
	share_cn_periodic(par->my_rank, Tij, in_r, in_l, par, tag);
	
	index = 0;
	for (i=0; i<=1-1; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_cn[start-1+i][j][k];
				Tij[ijk] = in_l[index];
				index += 1;
			}
		}
	}
	
	
	index = 0;
	for (i=0; i<=2-1; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_cn[end+i][j][k];
				Tij[ijk] = in_r[index];
				index += 1;
			}
		}
	}
	
	
	free_dvector(in_l, 0, 1*(ny+1)*nz-1);
	free_dvector(in_r, 0, 2*(ny+1)*nz-1);
	
	for (i=start-1; i<=end+1; ++i){
		for (k=0; k<=nz-1; ++k){
			
            // 4th order interpolation
			// T_0 = (- T_2 + 4*( T_1 + T_-1 ) - T_-2 )/6
			index = par->index_cn[i][-1][k];
			Tij[index]
			= 4.0*( Tij[par->index_cn[i][2][k]] + Tij[par->index_cn[i][0][k]] )
			- Tij[par->index_cn[i][3][k]]
			- 6.0*Tij[par->index_cn[i][1][k]];
            
			index = par->index_cn[i][ny+1][k];
			Tij[index]
			= 4.0*( Tij[par->index_cn[i][ny-2][k]] + Tij[par->index_cn[i][ny][k]] )
			- Tij[par->index_cn[i][ny-3][k]]
			- 6.0*Tij[par->index_cn[i][ny-1][k]];
			
			
		}
	}
	
	
	return;
}



void Tij_fill_cn__2(double *Tij, parameters *par, int tag)
{
	
	int i, j, k;
	int index;
	int ijk;
	const int ny = par->ny;
	const int nz = par->nz;
	const int end = par->local_nx_start + par->local_nx;
	
	double *in_r, *in_l;
	
	in_l = dvector(0, 1*(ny+1)*nz - 1); //freed
	in_r = dvector(0, 2*(ny+1)*nz - 1); //freed
	
	share_cn_periodic(par->my_rank, Tij, in_r, in_l, par, tag);
    
	index = 0;
	for (i=0; i<=1-1; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_cn[end+i][j][k];
				Tij[ijk] = in_r[index];
				index += 1;
			}
		}
	}
	
	
	free_dvector(in_l, 0, 1*(ny+1)*nz-1);
	free_dvector(in_r, 0, 2*(ny+1)*nz-1);
	
	return;
}

void Txx_fill_ct(double *Tij, parameters *par, int tag)
{
	int i, j, k;
	int index;
	int ijk;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *in_r, *in_l;
	
	in_l = dvector(0, 2*ny*nz - 1); //freed
	in_r = dvector(0, 2*ny*nz - 1); //freed
	
	share_ct_periodic(par->my_rank, Tij, in_r, in_l, par, tag);
	
	index = 0;
	for (i=0; i<=2-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_ct[start-2+i][j][k];
				Tij[ijk] = in_l[index];
				index += 1;
			}
		}
	}
	
	
	index = 0;
	for (i=0; i<=2-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_ct[end+i][j][k];
				Tij[ijk] = in_r[index];
				index += 1;
			}
		}
	}
	
	
	
	free_dvector(in_l, 0, 2*ny*nz-1);
	free_dvector(in_r, 0, 2*ny*nz-1);
	
	return;
}

void Txx_fill_ct__2(double *Tij, parameters *par, int tag)
{
	int i, j, k;
	int index;
	int ijk;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *in_r, *in_l;
	
	in_l = dvector(0, 2*ny*nz - 1); //freed
	in_r = dvector(0, 2*ny*nz - 1); //freed
	
	share_ct_periodic(par->my_rank, Tij, in_r, in_l, par, tag);
	
	index = 0;
	index = ny*nz;
	for (i=0; i<=1-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_ct[start-1+i][j][k];
				Tij[ijk] = in_l[index];
				index += 1;
			}
		}
	}
	
	
	index = 0;
	for (i=0; i<=1-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				ijk = par->index_ct[end+i][j][k];
				Tij[ijk] = in_r[index];
				index += 1;
			}
		}
	}
	
	
	
	free_dvector(in_l, 0, 2*ny*nz-1);
	free_dvector(in_r, 0, 2*ny*nz-1);
	
	return;
}


void Tyy_fill_ct(double *Tyy_cn, double *Tyy_ct, parameters *par)
{
	int i, k;
	int index;
	int ijk;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *bc_Tyy_t;
	double *bc_Tyy_b;
	double *Tyy_t1;
	double *Tyy_b1;
	double *scratch;
	
	
	scratch = FFTW_MALLOC(par->local_size_ved); //freed
	bc_Tyy_t= FFTW_MALLOC(local_nx*nz); //freed
	bc_Tyy_b= FFTW_MALLOC(local_nx*nz); //freed
	Tyy_t1= FFTW_MALLOC(local_nx*nz); //freed
	Tyy_b1= FFTW_MALLOC(local_nx*nz); //freed
	
	for (i=0; i<=par->local_size_ved-1; ++i){
		scratch[i] = 0.0;
	}
	
	interpolate_x1_cn2ved__4(Tyy_cn, scratch, par);
    
	index = 0;
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			ijk = par->index_ved[i][0][k];
			bc_Tyy_b[index] = 9.0/8.0*scratch[ijk];
			
			ijk = par->index_ved[i][ny][k];
			bc_Tyy_t[index] = 9.0/8.0*scratch[ijk];
			
			
			
			ijk = par->index_ved[i][1][k];
			Tyy_b1[index] = 9.0/8.0*scratch[ijk];
			
			ijk = par->index_ved[i][ny-1][k];
			Tyy_t1[index] = 9.0/8.0*scratch[ijk];
			
			
			index += 1;
		}
	}
    
	
	interpolate_x3_cn2ved__4(Tyy_cn, scratch, par);
	
	index = 0;
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			ijk = par->index_ved[i][0][k];
			bc_Tyy_b[index] += - 1.0/8.0*scratch[ijk];
			
			ijk = par->index_ved[i][ny][k];
			bc_Tyy_t[index] += - 1.0/8.0*scratch[ijk];
			
			
			ijk = par->index_ved[i][1][k];
			Tyy_b1[index] += - 1.0/8.0*scratch[ijk];
			
			ijk = par->index_ved[i][ny-1][k];
			Tyy_t1[index] += - 1.0/8.0*scratch[ijk];
			
			index += 1;
		}
	}
	
	
	index = 0;
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			ijk = par->index_ct[i][-1][k];
			Tyy_ct[ijk]
			= - 16.0*Tyy_b1[index]
			+ 9.0*( Tyy_ct[par->index_ct[i][1][k]] + Tyy_ct[par->index_ct[i][0][k]] )
			- 1.0*Tyy_ct[par->index_ct[i][2][k]];
			
			ijk = par->index_ct[i][ny][k];
			Tyy_ct[ijk]
			= - 16.0*Tyy_t1[index]
			+ 9.0*( Tyy_ct[par->index_ct[i][ny-1][k]] + Tyy_ct[par->index_ct[i][ny-2][k]] )
			- 1.0*Tyy_ct[par->index_ct[i][ny-3][k]];
			
			
			index += 1;
		}
	}
	
	
	fftw_free(scratch);
	fftw_free(bc_Tyy_b);
	fftw_free(bc_Tyy_t);
	fftw_free(Tyy_b1);
	fftw_free(Tyy_t1);
	
	return;	
}








