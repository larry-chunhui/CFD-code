/*
 * pentasolver.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Pentadiagonal matrix is defined (void get_pentamatrix_uw/v(double *D, int substep, parameters *par))
 * and the linear system is solved (void vel_pentasolver_u/v/w(double *u_scratch, double *rhs_u, int substep, parameters *par))
 * for the intermediate velocity (the first of the three fractional steps).
 * This is used for the fourth order cases.
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "pentasolver.h"



/* ******************************************************************* */

/* ******************************************************************* */
void vel_pentasolver_u(
                       double *u_scratch,
                       double *rhs_u,
                       int substep,
                       parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	
	double *D;
	double *scratch;
    
	int kl, ku, nrhs;
	int N, ldb, info;
    
	N = ny;
	nrhs = 1;
	kl = 2;
	ku = 2;
	int piv[N];
	ldb = 2*kl+ku+1;
	
	D = dvector(0, N*ldb-1);
	scratch = dvector(0, ny-1);
	
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
            
			// define right hand side
			for (j=0; j<=ny-1; ++j){
				scratch[j] = rhs_u[par->index_u[i][j][k]];
			}
			
			get_pentamatrix_uw(D, substep, par);
			
			dgbsv_(&N, &kl, &ku, &nrhs, D, &ldb, piv, scratch, &N, &info);
			
			for (j=0; j<=ny-1; ++j){
				u_scratch[par->index_u[i][j][k]] = scratch[j];
			}
			
		}
	}
	
	free_dvector(scratch, 0, ny-1);
	free_dvector(D, 0, N*ldb-1);
	
	return;
}


/* ******************************************************************* */

/* ******************************************************************* */
void vel_pentasolver_v(
                       double *v_scratch,
                       double *rhs_v,
                       int substep,
                       parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	
	double *D;
	double *scratch;
    
	int kl, ku, nrhs;
	int N, ldb, info;
    
	N = ny-1;
	nrhs = 1;
	kl = 2;
	ku = 2;
	int piv[N];
	ldb = 2*kl+ku+1;
	
	D = dvector(0, N*ldb-1);
	scratch = dvector(0, ny-2);
    
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
            
			// define right hand side
			for (j=1; j<=ny-1; ++j){
				scratch[j-1] = rhs_v[par->index_v[i][j][k]];
			}
			
			get_pentamatrix_v(D, substep, par);
			
			dgbsv_(&N, &kl, &ku, &nrhs, D, &ldb, piv, scratch, &N, &info);
			
			for (j=1; j<=ny-1; ++j){
				v_scratch[par->index_v[i][j][k]] = scratch[j-1];
			}
			
		}
	}
	
	free_dvector(scratch, 0, ny-2);
	free_dvector(D, 0, N*ldb-1);
	
	return;
}


/* ******************************************************************* */

/* ******************************************************************* */
void vel_pentasolver_w(
                       double *w_scratch,
                       double *rhs_w,
                       int substep,
                       parameters *par)
{
	
	int i, j, k;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	
	double *D;
	double *scratch;
    
	int kl, ku, nrhs;
	int N, ldb, info;
    
	N = ny;
	nrhs = 1;
	kl = 2;
	ku = 2;
	int piv[N];
	ldb = 2*kl+ku+1;
	
	D = dvector(0, N*ldb-1);
	scratch = dvector(0, N-1);
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
            
			// define right hand side
			for (j=0; j<=ny-1; ++j){
				scratch[j] = rhs_w[par->index_wp[i][j][k]];
			}
			
			get_pentamatrix_uw(D, substep, par);
			
			dgbsv_(&N, &kl, &ku, &nrhs, D, &ldb, piv, scratch, &N, &info);
			
			for (j=0; j<=ny-1; ++j){
				w_scratch[par->index_wp[i][j][k]] = scratch[j];
			}
			
		}
	}
	
	free_dvector(scratch, 0, ny-1);
	free_dvector(D, 0, N*ldb-1);
	
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void get_pentamatrix_uw(
                        double *D,
                        int substep,
                        parameters *par)
{
	
	int i, j, index;
	const int ny = par->ny;
	const double beta = par->rk3_beta[substep];
	const double coeff = par->dt*beta/par->re/( 12.0*pow(par->dy, 2.0) );
	
	const int kl = 2;
	const int ku = 2;
	const int ldb = 2*kl+ku+1; // 7
	
	index = 0;
	for (i=0; i<=ny-1; ++i){
		for (j=0; j<=ldb-1; ++j){
			D[index] = 0.0;
			index += 1;
		}
	}
	
	j = 0; // bottom wall
	D[j*(ldb)+2] = - 0.0*coeff;
	D[j*(ldb)+3] = - 0.0*coeff;
	D[j*(ldb)+4] = 1.0 - ( -53.0*coeff);
	D[j*(ldb)+5] = - ( 18.0*coeff);
	D[j*(ldb)+6] = - ( -1.0*coeff);
	
	j = 1;
	D[j*(ldb)+2] = - 0.0*coeff;
	D[j*(ldb)+3] = - (58.0/3.0)*coeff;
	D[j*(ldb)+4] = 1.0 - ( -(91.0/3.0)*coeff);
	D[j*(ldb)+5] = - ( 16.0*coeff );
	D[j*(ldb)+6] = - ( -1.0*coeff );
	
	
	for (j=2; j<=ny-3; ++j){
		D[j*(ldb)+2] = - ( -1.0*coeff );
		D[j*(ldb)+3] = - ( 16.0*coeff );
		D[j*(ldb)+4] = 1.0 - ( -30.0*coeff );
		D[j*(ldb)+5] = - ( 16.0*coeff );
		D[j*(ldb)+6] = - ( -1.0*coeff );
	}
	
	j = ny-2;
	D[j*(ldb)+6] = - 0.0*coeff;
	D[j*(ldb)+5] = - (58.0/3.0)*coeff;
	D[j*(ldb)+4] = 1.0 - ( -(91.0/3.0)*coeff);
	D[j*(ldb)+3] = - ( 16.0*coeff );
	D[j*(ldb)+2] = - ( -1.0*coeff );
	
	
	j = ny-1; // top wall
	D[j*(ldb)+6] = - 0.0*coeff;
	D[j*(ldb)+5] = - 0.0*coeff;
	D[j*(ldb)+4] = 1.0 - ( -53.0*coeff);
	D[j*(ldb)+3] = - ( 18.0*coeff);
	D[j*(ldb)+2] = - ( -1.0*coeff);
	
	
	return;
}



/* ******************************************************************* */

/* ******************************************************************* */
void get_pentamatrix_v(
                       double *D,
                       int substep,
                       parameters *par)
{
	
	int i, j, index;
	const int ny = par->ny;
	const int kl = 2;
	const int ku = 2;
	const int ldb = 2*kl+ku+1; // 7
	const double beta = par->rk3_beta[substep];
	const double coeff = par->dt*beta/par->re/( 12.0*pow(par->dy, 2.0) );
	
	index = 0;
	for (i=1; i<=ny-1; ++i){
		for (j=0; j<=ldb-1; ++j){
			D[index] = 0.0;
			index += 1;
		}
	}
	
	j = 1; // bottom wall
	D[(j-1)*(ldb)+2] = - ( 0.0*coeff );
	D[(j-1)*(ldb)+3] = - ( 0.0*coeff );
	D[(j-1)*(ldb)+4] = 1.0 - ( -29.0*coeff);
	D[(j-1)*(ldb)+5] = - ( 16.0*coeff);
	D[(j-1)*(ldb)+6] = - ( -1.0*coeff);
	
	j = 2;
	D[(j-1)*(ldb)+2] = - (  0.0*coeff );
	D[(j-1)*(ldb)+3] = - ( 16.0*coeff );
	D[(j-1)*(ldb)+4] = 1.0 - ( -30.0*coeff);
	D[(j-1)*(ldb)+5] = - ( 16.0*coeff);
	D[(j-1)*(ldb)+6] = - ( -1.0*coeff);
	
	
	for (j=3; j<=ny-3; ++j){
		D[(j-1)*(ldb)+2] = - ( -1.0*coeff );
		D[(j-1)*(ldb)+3] = - ( 16.0*coeff );
		D[(j-1)*(ldb)+4] = 1.0 - ( -30.0*coeff );
		D[(j-1)*(ldb)+5] = - ( 16.0*coeff );
		D[(j-1)*(ldb)+6] = - ( -1.0*coeff );
	}
	
	
	j = ny-2;
	D[(j-1)*(ldb)+6] = - (	0.0*coeff );
	D[(j-1)*(ldb)+5] = - ( 16.0*coeff );
	D[(j-1)*(ldb)+4] = 1.0 - ( -30.0*coeff);
	D[(j-1)*(ldb)+3] = - ( 16.0*coeff);
	D[(j-1)*(ldb)+2] = - ( -1.0*coeff);
	
	
	j = ny-1; // top wall
	D[(j-1)*(ldb)+6] = - ( 0.0*coeff );
	D[(j-1)*(ldb)+5] = - ( 0.0*coeff );
	D[(j-1)*(ldb)+4] = 1.0 - ( -29.0*coeff);
	D[(j-1)*(ldb)+3] = - ( 16.0*coeff );
	D[(j-1)*(ldb)+2] = - ( -1.0*coeff );
	
	
	return;
}
