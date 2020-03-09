/*
 *  share.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Exchange neighbouring data across different processes, including the periodic (streamwise) data.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "definitions.h"
#include "parameters.h"
#include "nrutil.h"
#include "share.h"


/* get left end data from my_rank = 0 to be used as the right boundary
 * at my_rank = nproc - 1
 * get right end data from my_rank = nproc - 1 to be used as the left boundary
 * at my_rank = 0
 */
void get_periodic_data(variables *var, neighbors *shared, parameters *par)
{
	int index, i, j, k;
	const int my_rank = par->my_rank;
	const int local_nx_start = par->local_nx_start;
	const int local_nx = par->local_nx;
	const int end = local_nx_start + local_nx;
	const int width = par->width_share;
	
	const int ny = par->ny;
	const int nz = par->nz;
	const int nproc = par->nproc;
	
	int tag;
	MPI_Status stat1, stat2, stat3;
	
	double *u_out;
	double *v_out;
	double *w_out;
    
	u_out = dvector(0, width*(ny+2)*nz - 1);
	v_out = dvector(0, width*(ny+1)*nz - 1);
	w_out = dvector(0, width*(ny+2)*nz - 1);
	
	if (my_rank == 0){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				u_out[index] = var->bc_u_b[par->index_tb[0+i][k]];
				w_out[index] = var->bc_w_b[par->index_tb[0+i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					u_out[index] = var->u[par->index_u[0+i][j][k]];
					w_out[index] = var->w[par->index_wp[0+i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				u_out[index] = var->bc_u_t[par->index_tb[0+i][k]];
				w_out[index] = var->bc_w_t[par->index_tb[0+i][k]];
				index = index + 1;
			}
		}
        
        
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = 0)
			for (k=0; k<=nz-1; ++k){
				v_out[index] = var->bc_v_b[par->index_tb[0+i][k]];
				index = index + 1;
			}
			// interior
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					v_out[index] = var->v[par->index_v[0+i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				v_out[index] = var->bc_v_t[par->index_tb[0+i][k]];
				index = index + 1;
			}
		}
		
		tag = 10;
		MPI_Sendrecv(u_out, width*(ny+2)*nz, MPI_DOUBLE, nproc-1, tag, \
                     shared->u_in_l, width*(ny+2)*nz, MPI_DOUBLE, nproc-1, tag, MPI_COMM_WORLD, &stat1);
		tag = 11;
		MPI_Sendrecv(v_out, width*(ny+1)*nz, MPI_DOUBLE, nproc-1, tag, \
                     shared->v_in_l, width*(ny+1)*nz, MPI_DOUBLE, nproc-1, tag, MPI_COMM_WORLD, &stat2);
		tag = 12;
		MPI_Sendrecv(w_out, width*(ny+2)*nz, MPI_DOUBLE, nproc-1, tag, \
                     shared->w_in_l, width*(ny+2)*nz, MPI_DOUBLE, nproc-1, tag, MPI_COMM_WORLD, &stat3);
        
	}
	
	
	if (my_rank == nproc - 1){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				u_out[index] = var->bc_u_b[par->index_tb[end - 1 - (width-1) + i][k]];
				w_out[index] = var->bc_w_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					u_out[index] = var->u[par->index_u[end - 1 - (width-1) + i][j][k]];
					w_out[index] = var->w[par->index_wp[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				u_out[index] = var->bc_u_t[par->index_tb[end - 1 - (width-1) + i][k]];
				w_out[index] = var->bc_w_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
        
        
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = 0)
			for (k=0; k<=nz-1; ++k){
				v_out[index] = var->bc_v_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					v_out[index] = var->v[par->index_v[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				v_out[index] = var->bc_v_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
		
		tag = 10;
		MPI_Sendrecv(u_out, width*(ny+2)*nz, MPI_DOUBLE, 0, tag, \
                     shared->u_in_r, width*(ny+2)*nz, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat1);
		tag = 11;
		MPI_Sendrecv(v_out, width*(ny+1)*nz, MPI_DOUBLE, 0, tag, \
                     shared->v_in_r, width*(ny+1)*nz, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat2);
		tag = 12;
		MPI_Sendrecv(w_out, width*(ny+2)*nz, MPI_DOUBLE, 0, tag, \
                     shared->w_in_r, width*(ny+2)*nz, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat3);
	}
	
	free_dvector(u_out, 0, width*(ny+2)*nz-1);
	free_dvector(v_out, 0, width*(ny+1)*nz-1);
	free_dvector(w_out, 0, width*(ny+2)*nz-1);
	
	return;
}



void get_neignboring_data(variables *var, neighbors *shared, parameters *par)
{
	const int ny = par->ny;
	const int nz = par->nz;
	const int width = par->width_share;
	
	shared->u_in_r = dvector(0, width*(ny+2)*nz-1);
	shared->u_in_l = dvector(0, width*(ny+2)*nz-1);
	shared->v_in_r = dvector(0, width*(ny+1)*nz-1);
	shared->v_in_l = dvector(0, width*(ny+1)*nz-1);
	shared->w_in_r = dvector(0, width*(ny+2)*nz-1);
	shared->w_in_l = dvector(0, width*(ny+2)*nz-1);
	
	
	share_u(par->my_rank, var, shared->u_in_r, shared->u_in_l, par, 0);
	share_v(par->my_rank, var, shared->v_in_r, shared->v_in_l, par, 1);
	share_w(par->my_rank, var, shared->w_in_r, shared->w_in_l, par, 2);
    
    
    
	return;
}

void finalize_neignboring(neighbors *shared, parameters *par)
{
    
	const int ny = par->ny;
	const int nz = par->nz;
	const int width = par->width_share;
	
	free_dvector(shared->u_in_r, 0, width*(ny+2)*nz-1);
	free_dvector(shared->u_in_l, 0, width*(ny+2)*nz-1);
	free_dvector(shared->v_in_r, 0, width*(ny+1)*nz-1);
	free_dvector(shared->v_in_l, 0, width*(ny+1)*nz-1);
	free_dvector(shared->w_in_r, 0, width*(ny+2)*nz-1);
	free_dvector(shared->w_in_l, 0, width*(ny+2)*nz-1);
    
	return;
}



void share_u(int my_rank, variables *var, double *u_in_r, double *u_in_l, parameters *par, int tag)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int width = par->width_share;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
	
	MPI_Status stat1, stat2;
    
	double *u_out_l;
	double *u_out_r;
    
	u_out_r = dvector(0, width*(ny+2)*nz - 1);
	u_out_l = dvector(0, width*(ny+2)*nz - 1);
	
	for (i=0; i<=width*(ny+2)*nz-1; ++i){
		u_out_r[i] = 0.0;
		u_out_l[i] = 0.0;
	}
	
	if (my_rank == 0){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boudary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				u_out_r[index] = var->bc_u_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					u_out_r[index] = var->u[par->index_u[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				u_out_r[index] = var->bc_u_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
		
		MPI_Sendrecv(u_out_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag,\
                     u_in_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat1);
		
		
	}else if(my_rank == par->nproc-1){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				u_out_l[index] = var->bc_u_b[par->index_tb[local_nx_start + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					u_out_l[index] = var->u[par->index_u[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				u_out_l[index] = var->bc_u_t[par->index_tb[local_nx_start + i][k]];
				index = index + 1;
			}
		}
		
		MPI_Sendrecv(u_out_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     u_in_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat2);
		
        
	}else{
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				u_out_l[index] = var->bc_u_b[par->index_tb[local_nx_start + i][k]];
				u_out_r[index] = var->bc_u_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					u_out_l[index] = var->u[par->index_u[local_nx_start + i][j][k]];
					u_out_r[index] = var->u[par->index_u[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				u_out_l[index] = var->bc_u_t[par->index_tb[local_nx_start + i][k]];
				u_out_r[index] = var->bc_u_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
		
        
		MPI_Sendrecv(u_out_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag,\
                     u_in_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		MPI_Sendrecv(u_out_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     u_in_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
	}
    
	
	free_dvector(u_out_r, 0, width*(ny+2)*nz-1);
	free_dvector(u_out_l, 0, width*(ny+2)*nz-1);
	
	return;
}





void share_v(int my_rank, variables *var, double *v_in_r, double *v_in_l, parameters *par, int tag)
{
	int i, j, k, index;
    const int ny = par->ny;
	const int nz = par->nz;
	const int width = par->width_share;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
    
	MPI_Status stat1, stat2;
    
	
	double *v_out_l;
	double *v_out_r;
	
	v_out_r = dvector(0, width*(ny+1)*nz - 1);
	v_out_l = dvector(0, width*(ny+1)*nz - 1);
	
	for (i=0; i<=width*(ny+1)*nz-1; ++i){
		v_out_r[i] = 0.0;
		v_out_l[i] = 0.0;
	}
	
	if (my_rank == 0){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = 0)
			for (k=0; k<=nz-1; ++k){
				v_out_r[index] = var->bc_v_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					v_out_r[index] = var->v[par->index_v[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				v_out_r[index] = var->bc_v_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
        
		MPI_Sendrecv(v_out_r, width*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     v_in_r, width*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat1);
		
		
	}else if(my_rank == par->nproc-1){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = 0)
			for (k=0; k<=nz-1; ++k){
				v_out_l[index] = var->bc_v_b[par->index_tb[local_nx_start + i][k]];
				index = index + 1;
			}
			// interior
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					v_out_l[index] = var->v[par->index_v[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				v_out_l[index] = var->bc_v_t[par->index_tb[local_nx_start + i][k]];
				index = index + 1;
			}
		}
        
		MPI_Sendrecv(v_out_l, width*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     v_in_l, width*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat2);
        
        
	}else{
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = 0)
			for (k=0; k<=nz-1; ++k){
				v_out_l[index] = var->bc_v_b[par->index_tb[local_nx_start + i][k]];
				v_out_r[index] = var->bc_v_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					v_out_l[index] = var->v[par->index_v[local_nx_start + i][j][k]];
					v_out_r[index] = var->v[par->index_v[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				v_out_l[index] = var->bc_v_t[par->index_tb[local_nx_start + i][k]];
				v_out_r[index] = var->bc_v_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
        
		MPI_Sendrecv(v_out_l, width*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     v_in_l, width*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		MPI_Sendrecv(v_out_r, width*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     v_in_r, width*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
        
	}
	
	
	free_dvector(v_out_r, 0, width*(ny+1)*nz-1);
	free_dvector(v_out_l, 0, width*(ny+1)*nz-1);
	
	return;
    
    
}



void share_w(int my_rank, variables *var, double *w_in_r, double *w_in_l, parameters *par, int tag)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int width = par->width_share;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
	
	MPI_Status stat1, stat2;
    
	
	double *w_out_l;
	double *w_out_r;
	
	w_out_r = dvector(0, width*(ny+2)*nz - 1);
	w_out_l = dvector(0, width*(ny+2)*nz - 1);
	
	
	for (i=0; i<=width*(ny+2)*nz-1; ++i){
		w_out_r[i] = 0.0;
		w_out_l[i] = 0.0;
	}
	
	if (my_rank == 0){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				w_out_r[index] = var->bc_w_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					w_out_r[index] = var->w[par->index_wp[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				w_out_r[index] = var->bc_w_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
        
		MPI_Sendrecv(w_out_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     w_in_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
        
		
	}else if(my_rank == par->nproc-1){
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				w_out_l[index] = var->bc_w_b[par->index_tb[local_nx_start + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					w_out_l[index] = var->w[par->index_wp[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				w_out_l[index] = var->bc_w_t[par->index_tb[local_nx_start + i][k]];
				index = index + 1;
			}
		}
        
		MPI_Sendrecv(w_out_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     w_in_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
        
	}else{
		index = 0;
		for (i=0; i<=width-1; ++i){
			// bottom boundary data (j = -1)
			for (k=0; k<=nz-1; ++k){
				w_out_l[index] = var->bc_w_b[par->index_tb[local_nx_start + i][k]];
				w_out_r[index] = var->bc_w_b[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
			// interior
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					w_out_l[index] = var->w[par->index_wp[local_nx_start + i][j][k]];
					w_out_r[index] = var->w[par->index_wp[end - 1 - (width-1) + i][j][k]];
					index = index + 1;
				}
			}
			// top boundary data (j = ny)
			for (k=0; k<=nz-1; ++k){
				w_out_l[index] = var->bc_w_t[par->index_tb[local_nx_start + i][k]];
				w_out_r[index] = var->bc_w_t[par->index_tb[end - 1 - (width-1) + i][k]];
				index = index + 1;
			}
		}
        
		MPI_Sendrecv(w_out_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     w_in_l, width*(ny+2)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		MPI_Sendrecv(w_out_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     w_in_r, width*(ny+2)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
        
	}
	
    
	
	free_dvector(w_out_r, 0, width*(ny+2)*nz-1);
	free_dvector(w_out_l, 0, width*(ny+2)*nz-1);
	
	return;
}




// index_cn[start-1 -> end+1][-1 -> ny+1][0 -> nz-1]
void share_cn(int my_rank, double *Tij, double *in_r, double *in_l, parameters *par, int tag)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
	
	MPI_Status stat1, stat2;
	
	
	double *out_l;
	double *out_r;
	
	out_r = dvector(0, 1*(ny+1)*nz - 1);
	out_l = dvector(0, 2*(ny+1)*nz - 1);
	
	
	for (i=0; i<=1*(ny+1)*nz-1; ++i){
		out_r[i] = 0.0;
	}
	for (i=0; i<=2*(ny+1)*nz-1; ++i){
		out_l[i] = 0.0;
	}
	
	
	if (my_rank == 0){
		index = 0;
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				out_r[index] = Tij[par->index_cn[end - 1][j][k]];
				index = index + 1;
			}
		}
		
		
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
					 in_r, 2*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}else if(my_rank == par->nproc-1){
		index = 0;
		for (i=0; i<=2-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_l[index] = Tij[par->index_cn[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_l, 2*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, \
					 in_l, 1*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		
	}else{
		index = 0;
		for (i=0; i<=2-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					out_l[index] = Tij[par->index_cn[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		index = 0;
		for (i=0; i<=1-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij[par->index_cn[end - 1][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_l, 2*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, \
					 in_l, 1*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
					 in_r, 2*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}
	
	
	
	free_dvector(out_r, 0, 1*(ny+1)*nz-1);
	free_dvector(out_l, 0, 2*(ny+1)*nz-1);
	
	return;
}





// index_cn[start-1 -> end+1][-1 -> ny+1][0 -> nz-1]
void share_cn_periodic(int my_rank, double *Tij, double *in_r, double *in_l, parameters *par, int tag)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
	
	MPI_Status stat1, stat2;
	
	
	double *out_l;
	double *out_r;
	
	out_r = dvector(0, 1*(ny+1)*nz - 1);
	out_l = dvector(0, 2*(ny+1)*nz - 1);
	
	
	for (i=0; i<=1*(ny+1)*nz-1; ++i){
		out_r[i] = 0.0;
	}
	for (i=0; i<=2*(ny+1)*nz-1; ++i){
		out_l[i] = 0.0;
	}
	
	
	if (my_rank == 0){
		index = 0;
		for (i=0; i<=1-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij[par->index_cn[end - 1][j][k]];
					index = index + 1;
				}
			}
		}
		
		index = 0;
		for (i=0; i<=2-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_l[index] = Tij[par->index_cn[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
        
		
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     in_r, 2*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
		MPI_Sendrecv(out_l, 2*(ny+1)*nz, MPI_DOUBLE, par->nproc-1, tag, \
                     in_l, 1*(ny+1)*nz, MPI_DOUBLE, par->nproc-1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}else if(my_rank == par->nproc-1){
		index = 0;
		for (i=0; i<=2-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_l[index] = Tij[par->index_cn[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		index = 0;
		for (i=0; i<=1-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij[par->index_cn[end-1][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_l, 2*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     in_l, 1*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, 0, tag, \
                     in_r, 2*(ny+1)*nz, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat1);
		
	}else{
		index = 0;
		for (i=0; i<=2-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_l[index] = Tij[par->index_cn[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		index = 0;
		for (i=0; i<=1-1; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij[par->index_cn[end - 1][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_l, 2*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, \
                     in_l, 1*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     in_r, 2*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}
	
	
	free_dvector(out_r, 0, 1*(ny+1)*nz-1);
	free_dvector(out_l, 0, 2*(ny+1)*nz-1);
	
	return;
}





void share_ct_periodic(int my_rank, double *Tij_ct, double *in_r, double *in_l, parameters *par, int tag)
{
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
	const int width = 2;
	
	MPI_Status stat1, stat2;
	
	
	double *out_l;
	double *out_r;
	
	out_r = dvector(0, width*ny*nz - 1);
	out_l = dvector(0, width*ny*nz - 1);
	
	
	for (i=0; i<=width*ny*nz-1; ++i){
		out_r[i] = 0.0;
		out_l[i] = 0.0;
	}
	
	if (my_rank == 0){
		index = 0;
		for (i=0; i<=width-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij_ct[par->index_ct[end-1 - (width-1) + i][j][k]];
					out_l[index] = Tij_ct[par->index_ct[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_r, width*ny*nz, MPI_DOUBLE, my_rank+1, tag, \
                     in_r, width*ny*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
		MPI_Sendrecv(out_l, width*ny*nz, MPI_DOUBLE, par->nproc-1, tag, \
                     in_l, width*ny*nz, MPI_DOUBLE, par->nproc-1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}else if(my_rank == par->nproc-1){
		index = 0;
		for (i=0; i<=width-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij_ct[par->index_ct[end-1 - (width-1) + i][j][k]];
					out_l[index] = Tij_ct[par->index_ct[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_l, width*ny*nz, MPI_DOUBLE, my_rank-1, tag, \
                     in_l, width*ny*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		
		MPI_Sendrecv(out_r, width*ny*nz, MPI_DOUBLE, 0, tag, \
                     in_r, width*ny*nz, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat1);
		
	}else{
		index = 0;
		for (i=0; i<=width-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					out_r[index] = Tij_ct[par->index_ct[end-1 - (width-1) + i][j][k]];
					out_l[index] = Tij_ct[par->index_ct[local_nx_start + i][j][k]];
					index = index + 1;
				}
			}
		}
		
		MPI_Sendrecv(out_l, width*ny*nz, MPI_DOUBLE, my_rank-1, tag, \
					 in_l, width*ny*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		MPI_Sendrecv(out_r, width*ny*nz, MPI_DOUBLE, my_rank+1, tag, \
					 in_r, width*ny*nz, MPI_DOUBLE, my_rank+1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}
	
	
	
	free_dvector(out_r, 0, width*ny*nz-1);
	free_dvector(out_l, 0, width*ny*nz-1);
	
	return;
}



void share_les_periodic__2(int my_rank, double *les, double *in_l, parameters *par, int tag)
{
	
	int i, j, k, index;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int end = local_nx_start + local_nx;
	
	MPI_Status stat1, stat2;
	
	double *out_r;
	
	out_r = dvector(0, 1*(ny+1)*nz - 1);
	
	for (i=0; i<=1*(ny+1)*nz-1; ++i){
		out_r[i] = 0.0;
		in_l[i]	 = 0.0;
	}
	
	
	index = 0;
	for (i=0; i<=1-1; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				out_r[index] = les[par->index_les[end-1][j][k]];
				index = index + 1;
			}
		}
	}
	
	if (my_rank == 0){
		
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     in_l, 1*(ny+1)*nz, MPI_DOUBLE, par->nproc-1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}else if(my_rank == par->nproc-1){
		
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, 0, tag, \
                     in_l, 1*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat1);
		
	}else{
		
		MPI_Sendrecv(out_r, 1*(ny+1)*nz, MPI_DOUBLE, my_rank+1, tag, \
                     in_l, 1*(ny+1)*nz, MPI_DOUBLE, my_rank-1, tag, MPI_COMM_WORLD, &stat2);
		
		
	}
	
	
	free_dvector(out_r, 0, 1*(ny+1)*nz-1);
	
	
	return;
	
}
