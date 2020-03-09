/*
 *  derivatives__4.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *  
 * Laplace, divergence and gradient calculations with fourth order schemes
 */

#include <mpi.h>
#include "definitions.h"
#include "parameters.h"
#include "nrutil.h"
#include "fft.h"
#include "derivatives__4.h"


/* ******************************************************************* */
/* boundary condition from Laplace operation in velocity               */
/* Lu_bc_y = vector(0, local_size_u-1)								   */
/* Lv_bc_y = vector(0, local_size_v-1)								   */
/* Lw_bc_y = vector(0, local_size_wp-1)								   */
/* ******************************************************************* */
void get_laplace_bc_y__4(variables *var, parameters *par)
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
	// interior ( excluding bottom and top cell)
	for (i=start; i<=end-1; ++i){
		for (j=2; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lu_y[par->index_u[i][j][k]] = 0.0;
			}
		}
	}
    
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 0;
			var->bc_Lu_y[par->index_u[i][j][k]] \
            = 104.0/3.0*var->bc_u_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = 1;
			var->bc_Lu_y[par->index_u[i][j][k]] \
            = - 8.0/3.0*var->bc_u_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-2;
			var->bc_Lu_y[par->index_u[i][j][k]] \
            = - 8.0/3.0*var->bc_u_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
            
			j = ny-1;
			var->bc_Lu_y[par->index_u[i][j][k]] \
            = 104.0/3.0*var->bc_u_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
		}
	}
	
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of v
	for (i=start; i<=end-1; ++i){
		for (j=3; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lv_y[par->index_v[i][j][k]] = 0.0;
			}
		}
	}
	
    
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 1;
			var->bc_Lv_y[par->index_v[i][j][k]] \
            = (16.0 - 2.0)*var->bc_v_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = 2;
			var->bc_Lv_y[par->index_v[i][j][k]] \
            = - 1.0*var->bc_v_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-2;
			var->bc_Lv_y[par->index_v[i][j][k]] \
            = - 1.0*var->bc_v_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-1;
			var->bc_Lv_y[par->index_v[i][j][k]] \
            = (16.0 - 2.0)*var->bc_v_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
		}
	}
	
    
	
	
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of w
	//interior
	for (i=start; i<=end-1; ++i){
		for (j=2; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lw_y[par->index_wp[i][j][k]] = 0.0;
			}
		}
	}
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 0;
			var->bc_Lw_y[par->index_wp[i][j][k]] \
            = 104.0/3.0*var->bc_w_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = 1;
			var->bc_Lw_y[par->index_wp[i][j][k]] \
            = - 8.0/3.0*var->bc_w_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-2;
			var->bc_Lw_y[par->index_wp[i][j][k]] \
            = - 8.0/3.0*var->bc_w_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
            
			j = ny-1;
			var->bc_Lw_y[par->index_wp[i][j][k]] \
            = 104.0/3.0*var->bc_w_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
		}
	}
	
	return;
}



/* ******************************************************************* */
/* laplace given u, v, w Gp, get Lu, Lv, Lw, LGp                       */
/* defined at the edge, same as velocity							   */
/* u = vector(0, local_size_u-1)									   */
/* Lu_y = vector(0, local_size_u-1)									   */
/* ******************************************************************* */
void get_laplace_y__4(double *u,
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
                = ( - 1.0*(u[par->index_ued[i][j-2][k]] + u[par->index_ued[i][j+2][k]]) \
                   +16.0*(u[par->index_ued[i][j-1][k]] + u[par->index_ued[i][j+1][k]]) \
                   -30.0* u[par->index_ued[i][j][k]])/(12.0*pow(dy, 2.0));
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
                = (	-     (v[par->index_ved[i][j-2][k]] + v[par->index_ved[i][j+2][k]]) \
                   +16.0*(v[par->index_ved[i][j-1][k]] + v[par->index_ved[i][j+1][k]]) \
                   -30.0* v[par->index_ved[i][j][k]])   /  (12.0*pow(dy, 2.0));
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
                = ( -     (w[par->index_ct[i][j-2][k]] + w[par->index_ct[i][j+2][k]]) \
                   +16.0*(w[par->index_ct[i][j-1][k]] + w[par->index_ct[i][j+1][k]]) \
                   -30.0* w[par->index_ct[i][j][k]])/(12.0*pow(dy, 2.0));
			}
		}
	}
	
	
	return;
}



/* ******************************************************************* */
/* boundary condition from Laplace operation in velocity               */
/* Lu_bc_xz = vector(0, nq_x-1)										   */
/* Lv_bc_xz = vector(0, nq_y-1)										   */
/* Lw_bc_xz = vector(0, nq_z-1)										   */
/* ******************************************************************* */
void get_laplace_bc_xz__4(variables *var, parameters *par)
{
	
	int i, j, k, index;
    const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
    
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of u
	if (my_rank == 0){
		
		// interior ( excluding bottom and top cell)
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					if (i == start){ //next to left boundary
						var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
                        
                        
					}else if(i == end-1){ //next to right boundary
						var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
						
                        
					}else{	//interior
						var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
						
					}
				}
			}
		}
		
		
        /* ******************************************************************* */
        // right end process
	}else if (my_rank == par->nproc - 1){
		
		// interior ( excluding bottom and top cell)
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					if (i == start){ //next to left boundary
						var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
                        
					}else if(i == end-1){ //next to right boundary
						var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
						
					}else{	//interior
						var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
                        
					}
				}
			}
		}
		
		
        
        
	}else{
        /* ******************************************************************* */
        // interior process
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					var->bc_Lu_xz[par->index_u[i][j][k]] = 0.0;
				}
			}
		}
		
	}
    
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of v
	if (my_rank == 0){
		
		for (i=start; i<=end-1; ++i){
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					if (i == start){ //left end cell
						var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;
                        
					}else if(i == end-1){ //right end cell
						var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;
                        
                        
					}else{	//interior
						var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;
						
						
					}
				}
			}
		}
        
        
		
		
        /* ******************************************************************* */
        // right end process
        // laplace of v
	}else if (my_rank == par->nproc - 1){
		
		for (i=start; i<=end-1; ++i){
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					if (i == start){ //left end cell
						var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;
                        
					}else if(i == end-1){ //right end cell
						var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;;
                        
						
					}else{	//interior
						var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;
                        
					}
				}
			}
		}
        
        
        
        /* ******************************************************************* */
        // interior process
        // laplace of v
	}else{
		
		for (i=start; i<=end-1; ++i){
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					var->bc_Lv_xz[par->index_v[i][j][k]] = 0.0;
				}
			}
		}
	}
	
    
	
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// left end process
	// laplace of w
	if (my_rank == 0){
		
		//interior
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					if ( i == start){
						var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
                        
					}else if ( i == end - 1){
						var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
                        
                        
					}else{
						var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
                        
					}
				}
			}
		}
        
        /* ******************************************************************* */
        // right end process
        // laplace of w
	}else if (my_rank == par->nproc - 1){
		
		//interior
		for (i=start; i<=end-1; ++i){
			index = nz;
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					if ( i == start){
						var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
                        
					}else if ( i == end - 1){
						var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
						
					}else{
						var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
                        
					}
				}
			}
		}
        
        
        /* ******************************************************************* */
        // interior process
        // laplace of w
	}else{
		
		//interior
		for (i=start; i<=end-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					var->bc_Lw_xz[par->index_wp[i][j][k]] = 0.0;
                    
				}
			}
		}
	}
    
    
	return;
}


/* ******************************************************************* */
/* laplace given u, v, w Gp, get Lu, Lv, Lw, LGp                       */
/* defined at the edge, same as velocity							   */
/* u = vector(0, nq_x-1)											   */
/* Lu_xz = vector(0, nq_x-1)										   */
/* ******************************************************************* */
void get_laplace_xz__4(double *u,
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
                = (- 1.0*(u[par->index_ued[i-2][j][k]] + u[par->index_ued[i+2][j][k]]) \
                   +16.0*(u[par->index_ued[i-1][j][k]] + u[par->index_ued[i+1][j][k]]) \
                   -30.0* u[par->index_ued[i][j][k]])/(12.0*pow(dx, 2.0)) \
                
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
                = ( -     (v[par->index_ved[i-2][j][k]] + v[par->index_ved[i+2][j][k]]) \
                   +16.0*(v[par->index_ved[i-1][j][k]] + v[par->index_ved[i+1][j][k]]) \
                   -30.0* v[par->index_ved[i][j][k]])  /   (12.0*pow(dx, 2.0)) \
                
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
                = ( -     (w[par->index_ct[i-2][j][k]] + w[par->index_ct[i+2][j][k]]) \
                   +16.0*(w[par->index_ct[i-1][j][k]] + w[par->index_ct[i+1][j][k]]) \
                   -30.0* w[par->index_ct[i][j][k]])/(12.0*pow(dx, 2.0)) \
                
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
/* boundary condition from Laplace operation in velocity               */
/* Lu_bc = vector(0, nq_x-1)										   */
/* Lv_bc = vector(0, nq_y-1)										   */
/* Lw_bc = vector(0, nq_z-1)										   */
/* ******************************************************************* */
void get_laplace_bc__4(variables *var, parameters *par)
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
	// interior ( excluding bottom and top cell)
	for (i=start; i<=end-1; ++i){
		for (j=2; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lu[par->index_u[i][j][k]] = 0.0;
			}
		}
	}
    
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 0;
			var->bc_Lu[par->index_u[i][j][k]] \
            = 104.0/3.0*var->bc_u_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = 1;
			var->bc_Lu[par->index_u[i][j][k]] \
            = - 8.0/3.0*var->bc_u_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-2;
			var->bc_Lu[par->index_u[i][j][k]] \
            = - 8.0/3.0*var->bc_u_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
            
			j = ny-1;
			var->bc_Lu[par->index_u[i][j][k]] \
            = 104.0/3.0*var->bc_u_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
		}
	}
	
    
    
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of v
	for (i=start; i<=end-1; ++i){
		for (j=3; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lv[par->index_v[i][j][k]] = 0.0;
			}
		}
	}
	
    
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 1;
			var->bc_Lv[par->index_v[i][j][k]] \
            = (16.0 - 2.0)*var->bc_v_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = 2;
			var->bc_Lv[par->index_v[i][j][k]] \
            = - 1.0*var->bc_v_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-2;
			var->bc_Lv[par->index_v[i][j][k]] \
            = - 1.0*var->bc_v_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-1;
			var->bc_Lv[par->index_v[i][j][k]] \
            = (16.0 - 2.0)*var->bc_v_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
		}
	}
	
    
	
	
	
	/* ******************************************************************* */
	/* ******************************************************************* */
	// laplace of w
	// interior
	for (i=start; i<=end-1; ++i){
		for (j=2; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Lw[par->index_wp[i][j][k]] = 0.0;
			}
		}
	}
	
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 0;
			var->bc_Lw[par->index_wp[i][j][k]] \
            = 104.0/3.0*var->bc_w_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = 1;
			var->bc_Lw[par->index_wp[i][j][k]] \
            = - 8.0/3.0*var->bc_w_b[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
			
			j = ny-2;
			var->bc_Lw[par->index_wp[i][j][k]] \
            = - 8.0/3.0*var->bc_w_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
            
			j = ny-1;
			var->bc_Lw[par->index_wp[i][j][k]] \
            = 104.0/3.0*var->bc_w_t[par->index_tb[i][k]]/(12.0*pow(dy, 2.0));
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
void get_laplace__4(double *u,
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
    const double dy = par->dy;
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
                = (- 1.0*(u[par->index_ued[i-2][j][k]] + u[par->index_ued[i+2][j][k]]) \
                   +16.0*(u[par->index_ued[i-1][j][k]] + u[par->index_ued[i+1][j][k]]) \
                   -30.0* u[par->index_ued[i][j][k]])/(12.0*pow(dx, 2.0)) \
                
                + ( - 1.0*(u[par->index_ued[i][j-2][k]] + u[par->index_ued[i][j+2][k]]) \
                   +16.0*(u[par->index_ued[i][j-1][k]] + u[par->index_ued[i][j+1][k]]) \
                   -30.0* u[par->index_ued[i][j][k]])/(12.0*pow(dy, 2.0)) \
                
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
                = ( -     (v[par->index_ved[i-2][j][k]] + v[par->index_ved[i+2][j][k]]) \
                   +16.0*(v[par->index_ved[i-1][j][k]] + v[par->index_ved[i+1][j][k]]) \
                   -30.0* v[par->index_ved[i][j][k]])  /   (12.0*pow(dx, 2.0)) \
                
                + (	-     (v[par->index_ved[i][j-2][k]] + v[par->index_ved[i][j+2][k]]) \
                   +16.0*(v[par->index_ved[i][j-1][k]] + v[par->index_ved[i][j+1][k]]) \
                   -30.0* v[par->index_ved[i][j][k]])   /  (12.0*pow(dy, 2.0)) \
                
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
                = ( -     (w[par->index_ct[i-2][j][k]] + w[par->index_ct[i+2][j][k]]) \
                   +16.0*(w[par->index_ct[i-1][j][k]] + w[par->index_ct[i+1][j][k]]) \
                   -30.0* w[par->index_ct[i][j][k]])/(12.0*pow(dx, 2.0)) \
                
                + ( -     (w[par->index_ct[i][j-2][k]] + w[par->index_ct[i][j+2][k]]) \
                   +16.0*(w[par->index_ct[i][j-1][k]] + w[par->index_ct[i][j+1][k]]) \
                   -30.0* w[par->index_ct[i][j][k]])/(12.0*pow(dy, 2.0)) \
                
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

/* ******************************************************************* */
void get_divergence_bc__4(variables *var,
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
	
    
	//interior
	for (i=start; i<=end-1; ++i){
		for (j=2; j<=ny-3; ++j){
			for (k=0; k<=nz-1; ++k){
				var->bc_Du[par->index_wp[i][j][k]] = 0.0;
			}
		}
	}
    
	// bottom boundary
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = 0;
			var->bc_Du[par->index_wp[i][j][k]] = -25.0*var->bc_v_b[par->index_tb[i][k]]/(24.0*dy);
			j = 1;
			var->bc_Du[par->index_wp[i][j][k]] =   1.0*var->bc_v_b[par->index_tb[i][k]]/(24.0*dy);
		}
	}
	
	// top boundary
	j = ny-1;
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			j = ny-1;
			var->bc_Du[par->index_wp[i][j][k]] = 25.0*var->bc_v_t[par->index_tb[i][k]]/(24.0*dy);
			j = ny-2;
			var->bc_Du[par->index_wp[i][j][k]] = -1.0*var->bc_v_t[par->index_tb[i][k]]/(24.0*dy);
		}
	}
    
    
	
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
void get_divergence__4(double *u,
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
	
	//interior
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				Du[par->index_wp[i][j][k]] \
                = (    1.0*(u[par->index_ued[i-1][j][k]] - u[par->index_ued[i+2][j][k]]) \
                   - 27.0*(u[par->index_ued[i][j][k]] - u[par->index_ued[i+1][j][k]]) ) \
                /(24.0*dx) \
                
                + (    1.0*(v[par->index_ved[i][j-1][k]] - v[par->index_ved[i][j+2][k]]) \
                   - 27.0*(v[par->index_ved[i][j][k]] - v[par->index_ved[i][j+1][k]]) ) \
                /(24.0*dy) \
                
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
void get_gradient__4(variables *var, int tag,
                     parameters *par)
{
    
	int i, j, k, index;
    const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
    const double dy = par->dy;
    const double dx = par->dx;
	
	
	int tag_s, tag_r;
	tag_r = par->my_rank;
	
	MPI_Status stat;
	
	double *pl_in;
	double *pl_out;
	double *pr_in;
	double *pr_out;
	pl_in = dvector(0, 2*ny*nz - 1);
	pl_out = dvector(0, 2*ny*nz - 1);
	pr_in = dvector(0, 2*ny*nz - 1);
	pr_out = dvector(0, 2*ny*nz - 1);
	
	for (i=0; i<=2*ny*nz-1; ++i){
		pl_in[i] = 0.0;
		pr_out[i] = 0.0;
	}
	
	for (i=0; i<=2*ny*nz-1; ++i){
		pr_in[i] = 0.0;
		pl_out[i] = 0.0;
	}
    
	
	index = 0;
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			pr_out[index] = var->p[par->index_wp[end - 2][j][k]];
			pl_out[index] = var->p[par->index_wp[start][j][k]];
			index = index + 1;
		}
	}
    
	for (j=0; j<=ny-1; ++j){
		for (k=0; k<=nz-1; ++k){
			pr_out[index] = var->p[par->index_wp[end - 1][j][k]];
			index = index + 1;
		}
	}
	
    
	if (my_rank == 0){
		
		tag_s = my_rank+1;
		MPI_Sendrecv(pr_out, 2*ny*nz, MPI_DOUBLE, my_rank+1, tag_s, \
					 pl_in, 2*ny*nz, MPI_DOUBLE, par->nproc-1, tag_r, MPI_COMM_WORLD, &stat);
		
		tag_s = par->nproc-1;
		MPI_Sendrecv(pl_out, ny*nz, MPI_DOUBLE, par->nproc-1, tag_s, \
					 pr_in, ny*nz, MPI_DOUBLE, my_rank+1, tag_r, MPI_COMM_WORLD, &stat);
        
	}else if(my_rank == par->nproc-1){
		
		tag_s = 0;
		MPI_Sendrecv(pr_out, 2*ny*nz, MPI_DOUBLE, 0, tag_s, \
					 pl_in, 2*ny*nz, MPI_DOUBLE, my_rank-1, tag_r, MPI_COMM_WORLD, &stat);
		
		tag_s = my_rank-1;
		MPI_Sendrecv(pl_out, ny*nz, MPI_DOUBLE, my_rank-1, tag_s, \
					 pr_in, ny*nz, MPI_DOUBLE, 0, tag_r, MPI_COMM_WORLD, &stat);
        
        
	}else{
		
		tag_s = my_rank+1;
		MPI_Sendrecv(pr_out, 2*ny*nz, MPI_DOUBLE, my_rank+1, tag_s, \
					 pl_in, 2*ny*nz, MPI_DOUBLE, my_rank-1, tag_r, MPI_COMM_WORLD, &stat);
		
		tag_s = my_rank-1;
		MPI_Sendrecv(pl_out, ny*nz, MPI_DOUBLE, my_rank-1, tag_s, \
					 pr_in, ny*nz, MPI_DOUBLE, my_rank+1, tag_r, MPI_COMM_WORLD, &stat);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
    
    
	// get dp/dx (= Gpx)
	// p defined at the cell center
	// dp/dx defined at the cell edge as velocity u
	// as dp/dx_(i, j) = (-p(i+1,j,k) + 27p(i,j,k) - 27p(i-1,j,k) + p(i-2,j,k))/dx
	for (i=start; i<=end-1; ++i){
		index = 0;
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (i == start){
					var->Gpx[par->index_u[i][j][k]] \
                    = ( pl_in[index] \
                       - 27.0*pl_in[index+ny*nz] \
                       + 27.0*var->p[par->index_wp[i][j][k]] \
                       - var->p[par->index_wp[i+1][j][k]] )/(24.0*dx);
                    
				}else if (i == start+1){
					var->Gpx[par->index_u[i][j][k]] \
                    = ( pl_in[index+ny*nz] \
                       - 27.0*var->p[par->index_wp[i-1][j][k]] \
                       + 27.0*var->p[par->index_wp[i][j][k]] \
                       - var->p[par->index_wp[i+1][j][k]] )/(24.0*dx);
                    
				}else if (i == end-1){
					var->Gpx[par->index_u[i][j][k]] \
                    = ( var->p[par->index_wp[i-2][j][k]] \
                       - 27.0*var->p[par->index_wp[i-1][j][k]] \
                       + 27.0*var->p[par->index_wp[i][j][k]] \
                       - pr_in[index] )/(24.0*dx);
				}else{
					var->Gpx[par->index_u[i][j][k]] \
                    = ( var->p[par->index_wp[i-2][j][k]] \
                       - 27.0*var->p[par->index_wp[i-1][j][k]] \
                       + 27.0*var->p[par->index_wp[i][j][k]] \
                       - var->p[par->index_wp[i+1][j][k]] )/(24.0*dx);
				}
				
				index += 1;
			}
		}
	}
	
	
	
	// get dp/dy (=Gpy)
	// p defined at the cell center
	// dp/dy defined at the cell edge as velocity v
	// as dp/dy_(i, j) = (-p(i,j+1,k) + 27p(i,j,k) - 27p(i,j-1,k) + p(i,j-2,k))/24dy
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (j == 1){
					var->Gpy[par->index_v[i][j][k]] \
                    = ( - var->p[par->index_wp[i][j+1][k]] \
                       + 27.0*var->p[par->index_wp[i][j][k]] - 26.0*var->p[par->index_wp[i][j-1][k]] )/(24.0*dy);
				}else if (j == ny-1){
					var->Gpy[par->index_v[i][j][k]] \
                    = ( var->p[par->index_wp[i][j-2][k]]  \
                       - 27.0*var->p[par->index_wp[i][j-1][k]] + 26.0*var->p[par->index_wp[i][j][k]] )/(24.0*dy);
                }else{
					var->Gpy[par->index_v[i][j][k]] \
                    = ( (var->p[par->index_wp[i][j-2][k]] - var->p[par->index_wp[i][j+1][k]]) \
                       - 27.0*(var->p[par->index_wp[i][j-1][k]] - var->p[par->index_wp[i][j][k]]) )/(24.0*dy);
				}
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
	
	
	free_dvector(pl_in, 0, 2*ny*nz - 1);
	free_dvector(pr_out, 0, 2*ny*nz - 1);
	free_dvector(pr_in, 0, 2*ny*nz - 1);
	free_dvector(pl_out, 0, 2*ny*nz - 1);
    
	return;
}
