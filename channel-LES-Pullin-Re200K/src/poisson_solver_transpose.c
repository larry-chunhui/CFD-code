/*
 *  poisson_solver.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * The second of the three fractional steps -- solving for the pressure.
 * Data is already in the frequency space in z. In addition, Cosine transform in y and Fourier transform in x.
 * Details on the numerical methods is found in 
 * M Inoue. Large-Eddy Simulation of the Flat-plate Turbulent Boundary
 * Layer at High Reynolds numbers. PhD thesis, California Institute of Technology, 2012.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "poisson_solver_transpose.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

#define POISSON_y(a,b,c) (b-0) + par->ny*(c-0) + par->ny*par->nz*(a-par->local_nx_start)
// Transpose method
#define TRANS_POISSON(a,b,c) (a-0) + par->nx*(b-0) + par->nx*par->ny*(c-start_z)


/* ******************************************************************* */
/* Solving pressure equation using cosine transform in y-direction	   */
/* i.e. solving helmholtz equation for pressure						   */
/* ******************************************************************* */
void poisson_1DDCT_in_y_transpose_1DFFT_in_x(
                                             double *delta_p,
                                             double *rhs_b,
                                             parameters *par,
                                             fftwplans *planptr)
{
	
	
	int i, j, k, index;
	const int my_rank = par->my_rank;
	const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
    const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	const int local_nz = nz/par->nproc;
	const int start_z = local_nz*my_rank;
	const int end_z = local_nz*(my_rank+1);
	const int local_size_wp_yzplane = par->local_size_wp;
	const int local_size_wp_xyplane = ny*nx*local_nz;
	
	
	if (nz%par->nproc != 0){
		printf("poisson_1DDCT_in_y_transpose(): nz has to be a multiple of par->nproc\n");
	}
	
	double *scratch_wp;
	double *rhs_bhat;
	double *delta_phat;
    
	scratch_wp = FFTW_MALLOC(local_size_wp_yzplane);
	delta_phat = FFTW_MALLOC(local_size_wp_xyplane);
	rhs_bhat = FFTW_MALLOC(local_size_wp_xyplane);
	
	double *scratch_xy;
    double *scratch_xy2;
	scratch_xy = FFTW_MALLOC(local_size_wp_xyplane);
    scratch_xy2 = FFTW_MALLOC(local_size_wp_xyplane);
	for (i=0; i<=local_size_wp_xyplane-1; ++i){
		scratch_xy[i] = 0.0; // initialize
        scratch_xy2[i] = 0.0; // initialize
	}
    
	// arrage a data structure in [i][k][j] format
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				scratch_wp[POISSON_y(i,j,k)] = rhs_b[par->index_wp[i][j][k]];
			}
		}
	}
	
	
    
	// Cosine transform in y
    fftw_forward_poisson_ver2(scratch_wp, scratch_wp, planptr, par);
	
	
	
	// transpose data in [k][i][j] format
	fftw_execute_r2r(planptr->transpose_yz2xy, scratch_wp, scratch_xy);
	
	
	// arrange a data structure in [k][j][i] format
	index = 0;
	for (k=start_z; k<=end_z-1; ++k){
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				scratch_xy2[TRANS_POISSON(i,j,k)] = scratch_xy[index];
				index += 1;
			}
		}
	}
    
    // Fourier transform in x
    fftw_execute_r2r(planptr->p1d_x_wp, scratch_xy2, rhs_bhat);
    double invnx = 1.0/(double)nx;
    for ( i=0; i<=local_size_wp_xyplane-1; ++i){
        rhs_bhat[i] = rhs_bhat[i]*invnx;
    }
	
	if (par->fd_order == 2){
		trisolver_helmholz_transpose_fft(scratch_xy2, rhs_bhat, par);
	}
	if (par->fd_order == 4){
		septasolver_helmholz_transpose_fft(scratch_xy2, rhs_bhat, par);
	}
    
	
    // Inverse Fourier transform in x
    fftw_execute_r2r(planptr->p1d_invx_wp, scratch_xy2, delta_phat);
	
	
	// Rearrange a data structure of the solution
    // in [k][i][j] format
	index = 0;
	for (k=start_z; k<=end_z-1; ++k){
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				scratch_xy[index] = delta_phat[TRANS_POISSON(i,j,k)];
				index += 1;
			}
		}
	}
	
    // transpose data in  [i][k][j] format
	fftw_execute_r2r(planptr->transpose_xy2yz, scratch_xy, scratch_wp);
	
	
	// Cosine inverse transform in y
	fftw_inverse_poisson_ver2(scratch_wp, scratch_wp, planptr, par);
	
	
	// Rearrange a data structure to finally come back to its original form
	// thread division in x
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				delta_p[par->index_wp[i][j][k]] = scratch_wp[POISSON_y(i,j,k)];
			}
		}
	}
	
	
	
	fftw_free(delta_phat);
	fftw_free(rhs_bhat);
	fftw_free(scratch_wp);
	fftw_free(scratch_xy);
    fftw_free(scratch_xy2);
	
	return;
}


/* ******************************************************************* */
/* Tridiagonal solver for helmholz equation of pressure				   */
/* Using lapack, pthrea-only environment							   */
/* ******************************************************************* */
void trisolver_helmholz_transpose_fft(
                                      double *delta_phat,
                                      double *rhs_bhat,
                                      parameters *par)
{
	int i, j, k;
	const int my_rank = par->my_rank;
	const int nx = par->nx;
	const int local_nz = par->nz/par->nproc;
	const int ny = par->ny;
    const int nz = par->nz;
	const int start_z = local_nz*my_rank;
	const int end_z = local_nz*(my_rank+1);
    const double dy = par->dy;
    const double dy2 = dy*dy;
    const double dx = par->dx;
    const double dx2 = dx*dx;
    
    
    for (i=0; i<=nx-1; ++i){
        for (j=0; j<=ny-1; ++j){
            for (k=start_z; k<=end_z-1; ++k){
                
                const double phase_x = 2*PI*i/(double)nx;
                const double element_x = 2*(cos(phase_x) - 1)/dx2;
                
                const double phase_y = PI*j/(double)ny;
                const double element_y = 2*(cos(phase_y) - 1)/dy2;
                const double element_z = - pow(par->kz[k], 2.0);
                
                double hoge;
                // singular equation
                if (i == 0 && (k == 0 || k == nz/2) && j == 0){
                    hoge = 0.0;
                }else{
                    hoge = rhs_bhat[TRANS_POISSON(i,j,k)];
                    hoge *= 1.0/(element_x+element_y+element_z);
                }
                
                delta_phat[TRANS_POISSON(i,j,k)] = hoge;
            }
        }
    }
	
	
	return;
}



/* ******************************************************************* */
/* Septadiagonal solver for helmholz equation of pressure              */
/* Using lapack, pthrea-only environment							   */
/* ******************************************************************* */
void septasolver_helmholz_transpose_fft(
                                        double *delta_phat,
                                        double *rhs_bhat,
                                        parameters *par)
{
    int i, j, k;
	const int my_rank = par->my_rank;
	const int nx = par->nx;
	const int local_nz = par->nz/par->nproc;
	const int ny = par->ny;
    const int nz = par->nz;
	const int start_z = local_nz*my_rank;
	const int end_z = local_nz*(my_rank+1);
    const double dy = par->dy;
    const double dy2 = dy*dy;
    const double dx = par->dx;
    const double dx2 = dx*dx;
    
    for (i=0; i<=nx-1; ++i){
        for (j=0; j<=ny-1; ++j){
            for (k=start_z; k<=end_z-1; ++k){
                
                const double phase_x = 2*PI*i/(double)nx;
                const double element_x = 2*( cos(3*phase_x) - 54*cos(2*phase_x)
                                            + 783*cos(phase_x) - 730)/(576*dx2);
                
                const double phase_y = PI*j/(double)ny;
                const double element_y = 2*( cos(3*phase_y) - 54*cos(2*phase_y)
                                            + 783*cos(phase_y) - 730)/(576*dy2);
                const double element_z = - pow(par->kz[k], 2.0);
                
                double hoge;
                // singular equation
                if (i == 0 && (k == 0 || k == nz/2) && j == 0){
                    hoge = 0.0;
                }else{
                    hoge = rhs_bhat[TRANS_POISSON(i,j,k)];
                    hoge *= 1.0/(element_x+element_y+element_z); 
                }
                
                delta_phat[TRANS_POISSON(i,j,k)] = hoge;
            }
        }
    }
	
	
	return;
}
