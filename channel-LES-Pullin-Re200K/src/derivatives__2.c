/*
 *  derivatives__2.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Laplace, divergence and gradient calculations with second order schemes.
 */

#include <mpi.h>
#include "definitions.h"
#include "parameters.h"
#include "nrutil.h"
#include "fft.h"
#include "derivatives__2.h"


/* ******************************************************************* */
/* boundary condition from Laplace operation in velocity               */
/* Lu_bc_y = vector(0, nq_x-1)										   */
/* Lv_bc_y = vector(0, nq_y-1)										   */
/* Lw_bc_y = vector(0, nq_z-1)										   */
/* ******************************************************************* */
void get_laplace_bc_y__2(variables *var, parameters *par)
{
	
	int i, j, k, index;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
    const double dy = par->dy;
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of u
	
	// interior ( excluding bottom and top cell)
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-2; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lu_y[par->index_u[i][j][k]] = 0.0;
			}
		}
	}
    
	j = 0;  //bottom end cell
	for (i=start; i<=end-1; ++i){
		index = nz*j;
		for (k=0; k<=nz-1; ++k){
            var->bc_Lu_y[par->index_u[i][j][k]] \
            = 2.0*var->bc_u_b[par->index_tb[i][k]]/pow(dy, 2.0);
		}
	}
	
	j = ny-1; //top end cell
	for (i=start; i<=end-1; ++i){
		index = nz*j;
		for (k=0; k<=nz-1; ++k){
            var->bc_Lu_y[par->index_u[i][j][k]] \
            = 2.0*var->bc_u_t[par->index_tb[i][k]]/pow(dy, 2.0);
		}
	}
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of v
	
	for (i=start; i<=end-1; ++i){
		index = nz;
		for (j=2; j<=ny-2; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lv_y[par->index_v[i][j][k]] = 0.0;
			}
		}
	}
	
	
	j = 1; //next to bottom boundary
	for (i=start; i<=end-1; ++i){
		index = nz*(j-1);
		for (k=0; k<=nz-1; ++k){
			var->bc_Lv_y[par->index_v[i][j][k]] \
            = var->bc_v_b[par->index_tb[i][k]]/pow(dy, 2.0);
		}
	}
    
	
	j = ny-1; //next to top boundary
	for (i=start; i<=end-1; ++i){
		index = nz*(j-1);
		for (k=0; k<=nz-1; ++k){
			var->bc_Lv_y[par->index_v[i][j][k]] \
            = var->bc_v_t[par->index_tb[i][k]]/pow(dy, 2.0);
		}
	}
	
	
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of w
	//interior
	for (i=start; i<=end-1; ++i){
		index = nz;
		for (j=1; j<=ny-2; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lw_y[par->index_wp[i][j][k]] = 0.0;
			}
		}
	}
	
	
	
	// bottom end cell except for left and right end
	j = 0;
	for (i=start; i<=end-1; ++i){
		index = nz*j;
		for (k=0; k<=nz-1; ++k){
			var->bc_Lw_y[par->index_wp[i][j][k]] \
            = 2.0*var->bc_w_b[par->index_tb[i][k]]/pow(dy, 2.0);
		}
	}
	
	
	
	// top end cell except for left and right end
	j = ny-1;
	for (i=start; i<=end-1; ++i){
		index = nz*j;
		for (k=0; k<=nz-1; ++k){
			var->bc_Lw_y[par->index_wp[i][j][k]] \
            = 2.0*var->bc_w_t[par->index_tb[i][k]]/pow(dy, 2.0);
		}
	}
	
    
	return;
}



/* ******************************************************************* */
/* laplace given u, v, w Gp, get Lu, Lv, Lw, LGp                       */
/* defined at the edge, same as velocity							   */
/* u = vector(0, nq_x-1)											   */
/* Lu_y = vector(0, nq_x-1)											   */
/* ******************************************************************* */
void get_laplace_y__2(double *u,
                      double *Lu,
                      double *v,
                      double *Lv,
                      double *w,
                      double *Lw,
                      variables *var,
                      parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
    const double dy = par->dy;
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of u
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Lu[par->index_u[i][j][k]] \
                = (u[par->index_ued[i][j+1][k]]- \
				   2.0*u[par->index_ued[i][j][k]]+ \
				   u[par->index_ued[i][j-1][k]])/pow(dy, 2.0);
			}
		}
	}
	
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of v
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Lv[par->index_v[i][j][k]] = (v[par->index_ved[i][j+1][k]]- \
											 2.0*v[par->index_ved[i][j][k]]+ \
											 v[par->index_ved[i][j-1][k]])/pow(dy, 2.0);
			}
		}
	}
	
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of w
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Lw[par->index_wp[i][j][k]] \
                = (w[par->index_ct[i][j+1][k]]-2.0*w[par->index_ct[i][j][k]]+w[par->index_ct[i][j-1][k]])/pow(dy, 2.0);
			}
		}
	}
	
	
	return;
}


/* ******************************************************************* */
/* laplace given u, v, w Gp, get Lu, Lv, Lw, LGp                       */
/* defined at the edge, same as velocity							   */
/* u = vector(0, nq_x-1)											   */
/* Lu = vector(0, nq_x-1)											   */
/* ******************************************************************* */
void get_laplace_xz__2(double *u,
                       double *Lu,
                       double *v,
                       double *Lv,
                       double *w,
                       double *Lw,
                       variables *var,
                       parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
    const int local_size_u = par->local_size_u;
    const int local_size_v = par->local_size_v;
    const int local_size_wp = par->local_size_wp;
    const double dx = par->dx;
	
	
	double *d2udz2;
	d2udz2 = FFTW_MALLOC(local_size_u);
	for (i=0; i<=local_size_u-1; ++i){
		d2udz2[i] = 0.0;
	}
	
	double *d2vdz2;
	d2vdz2 = FFTW_MALLOC(local_size_v);
	for (i=0; i<=local_size_v-1; ++i){
		d2vdz2[i] = 0.0;
	}
	
	double *d2wdz2;
	d2wdz2 = FFTW_MALLOC(local_size_wp);
	for (i=0; i<=local_size_wp-1; ++i){
		d2wdz2[i] = 0.0;
	}
    
    
    
	// d2u/dz2
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				d2udz2[par->index_u[i][j][k]] =  - var->u[par->index_u[i][j][k]]*pow(par->kz[k], 2.0);
			}
            
		}
	}
	
	// d2v/d2z
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				d2vdz2[par->index_v[i][j][k]] =  - var->v[par->index_v[i][j][k]]*pow(par->kz[k], 2.0);
			}
            
		}
	}
	
	// d2w/d2z
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				d2wdz2[par->index_wp[i][j][k]] =  - var->w[par->index_wp[i][j][k]]*pow(par->kz[k], 2.0);
			}
            
		}
	}
    
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of u
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Lu[par->index_u[i][j][k]] \
                = (u[par->index_ued[i+1][j][k]]-2.0*u[par->index_ued[i][j][k]]+u[par->index_ued[i-1][j][k]])/pow(dx, 2.0) \
                + d2udz2[par->index_u[i][j][k]];
			}
		}
	}
	
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of v
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Lv[par->index_v[i][j][k]] \
                = (v[par->index_ved[i+1][j][k]]-2.0*v[par->index_ved[i][j][k]]+v[par->index_ved[i-1][j][k]])/pow(dx, 2.0) \
                + d2vdz2[par->index_v[i][j][k]];
			}
		}
	}
	
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of w
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Lw[par->index_wp[i][j][k]] \
                = (w[par->index_ct[i+1][j][k]]-2.0*w[par->index_ct[i][j][k]]+w[par->index_ct[i-1][j][k]])/pow(dx, 2.0) \
                + d2wdz2[par->index_wp[i][j][k]];
			}
		}
	}
	
	
	fftw_free(d2udz2);
	fftw_free(d2vdz2);
	fftw_free(d2wdz2);
	
	return;
}



/* ******************************************************************* */
/* divergence of u or Gp : given u (or Gp) get Du or DGp			   */
/* u (Gpx) = vector(0, nq_x-1)	defined at the edge					   */
/* v (Gpy) = vector(0, nq_y-1)	defined at the edge					   */
/* w (Gpz) = vector(0, nq_z-1)	defined at the center				   */
/* Du = vector(0, np-1)												   */
/* DGp = vector(0, np-1) defined at the center  					   */
/* ******************************************************************* */
void get_divergence__2(double *u,
                       double *v,
                       double *w,
                       double *Du,
                       parameters *par)
{
	
	int i, j, k;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
    const int local_size_wp = par->local_size_wp;
    const double dy = par->dy;
    const double dx = par->dx;
	
	
	double *dwdz;
	dwdz = FFTW_MALLOC(local_size_wp);
	
	for (i=0; i<=local_size_wp-1; ++i){
		dwdz[i] = 0.0;
	}
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=(nz-1)/2; ++k){
				if (k == 0){
					dwdz[par->index_wp[i][j][k]] = 0.0;
				}else{
					dwdz[par->index_wp[i][j][nz-k]] = w[par->index_wp[i][j][k]]*par->kz[k];
					dwdz[par->index_wp[i][j][k]] =  - w[par->index_wp[i][j][nz-k]]*par->kz[nz-k];
				}
				
			}
		}
	}
    
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Du[par->index_wp[i][j][k]] \
                = (u[par->index_ued[i+1][j][k]] - u[par->index_ued[i][j][k]])/dx \
                + (v[par->index_ved[i][j+1][k]] - v[par->index_ved[i][j][k]])/dy \
                + dwdz[par->index_wp[i][j][k]];
			}
		}
	}
	
    
	
	fftw_free(dwdz);
    
	return;
}



/* ******************************************************************* */
/* gradient of p : given p get Gpx, Gpy and Gpz 					   */
/* p = vector(0, np-1)											       */
/* Gpx = vector(0, nq_x-1)											   */
/* Gpy = vector(0, nq_y-1)											   */
/* Gpz = vector(0, nq_z-1)											   */
/* ******************************************************************* */
void get_gradient__2(variables *var, int tag,
                     parameters *par)
{
    
	int i, j, k, index;
	int my_rank; my_rank = par->my_rank;
	int local_nx; local_nx = par->local_nx;
	int local_nx_start; local_nx_start = par->local_nx_start;
	int ny; ny = par->ny;
	int nz; nz = par->nz;
	int local_size_u; local_size_u = par->local_size_u;
	int local_size_v; local_size_v = par->local_size_v;
	int local_size_wp; local_size_wp = par->local_size_wp;
	int start; start = local_nx_start;
	int end; end = local_nx_start + local_nx;
	
	int mz; mz = par->mz;
	
	double dx; dx = par->dx;
	double dy; dy = par->dy;
	
	MPI_Status stat;
	
	double *p_in;
	double *p_out;
	p_in = dvector(0, ny*nz - 1);
	p_out = dvector(0, ny*nz - 1);
	
	for (i=0; i<=ny*nz-1; ++i){
		p_in[i] = 0.0;
		p_out[i] = 0.0;
	}
	
	if ( par->nproc != 1){
		
		if (my_rank == 0){
			index = 0;
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					p_out[index] = var->p[par->index_wp[end - 1][j][k]];
					index = index + 1;
				}
			}
			MPI_Sendrecv(p_out, ny*nz, MPI_DOUBLE, my_rank+1, tag, p_in, ny*nz, MPI_DOUBLE, par->nproc-1, tag, MPI_COMM_WORLD, &stat);
            
		}else if(my_rank == par->nproc-1){
			index = 0;
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					p_out[index] = var->p[par->index_wp[end - 1][j][k]];
					index = index + 1;
				}
			}
			MPI_Sendrecv(p_out, ny*nz, MPI_DOUBLE, 0, tag, p_in, ny*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat);
            
		}else{
			index = 0;
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					p_out[index] = var->p[par->index_wp[end - 1][j][k]];
					index = index + 1;
				}
			}
			MPI_Sendrecv(p_out, ny*nz, MPI_DOUBLE, my_rank+1, tag, p_in, ny*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat);
			
		}
	}
    
	
	// get dp/dx (= Gpx)
	// p defined at the cell center
	// dp/dx defined at the cell edge as velocity u
	// as dp/dx_(i, j) = (p(i,j,k) - p(i-1,j,k))/dx
	for (i=start; i<=end-1; ++i){
		index = 0;
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (i == start){
					var->Gpx[par->index_u[i][j][k]] = (var->p[par->index_wp[i][j][k]] - p_in[index])/dx;
					index = index + 1;
				}else{
					var->Gpx[par->index_u[i][j][k]] = (var->p[par->index_wp[i][j][k]] - var->p[par->index_wp[i-1][j][k]])/dx;
				}
			}
		}
	}
	
	
	// get dp/dy (=Gpy)
	// p defined at the cell center
	// dp/dy defined at the cell edge as velocity v
	// as dp/dy_(i, j) = (p(i,j,k) - p(i,j-1,k))/dy
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				var->Gpy[par->index_v[i][j][k]] = (var->p[par->index_wp[i][j][k]] - var->p[par->index_wp[i][j-1][k]])/dy;
			}
		}
	}
	
	// get dp/dz (=Gpz)
	// p defined at the cell center
	// dp/dz defined at the cell center as velocity w
	// as dp/dz_(i, j) = -i*kz*p(i,j,k)
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=(nz-1)/2; ++k){
				if (k == 0){
					var->Gpz[par->index_wp[i][j][k]] = 0.0;
				}else{
					var->Gpz[par->index_wp[i][j][k]] = - var->p[par->index_wp[i][j][nz-k]]*par->kz[nz-k];
					var->Gpz[par->index_wp[i][j][nz-k]] = var->p[par->index_wp[i][j][k]]*par->kz[k];
				}
				
			}
		}
	}
    
	free_dvector(p_in, 0, ny*nz - 1);
	free_dvector(p_out, 0, ny*nz - 1);
    
	return;
}
