/*
 *  parameters.c
 *
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Parameter reading, defining, allocating and deallocating.
 * Indices for different variables are also defined here.
 */


#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"


static void readparam(parameters *par);
extern void init_par(
                     parameters *par)
{
	
	MPI_Comm_rank(MPI_COMM_WORLD, &par->my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &par->nproc);
	
	readparam(par);
	
	int i, j, k, index;
	const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nx_start = par->local_nx_start;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	const int width = par->width_share;
	int width_v;
	int width_uw;
	
	if (par->fd_order == 2){
		width_v = 1;
		width_uw = 1;
	}
	if (par->fd_order == 4){
		width_v = 3;
		width_uw = 2;
	}
	
	par->index_u = i3tensor(start, end-1, 0, ny-1, 0, nz-1); // freed
	par->index_v = i3tensor(start, end-1, 1, ny-1, 0, nz-1); // freed
	par->index_wp = i3tensor(start, end-1, 0, ny-1, 0, nz-1); // freed
	
	if (par->my_rank == 0){
		
		index = 0;
		par->index_bigwp = i3tensor(0, nx-1, 0, ny-1, 0, nz-1);
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					par->index_bigwp[i][j][k] = index;
					index += 1;
				}
			}
		}
        
	}
	
	
	if (par->fd_order == 2){
		par->index_cn = i3tensor(start, end, 0, ny, 0, nz-1);
		par->index_ct = i3tensor(start-width, end-1+width, -1, ny, 0, nz-1);
		par->index_ued = i3tensor(start-width, end-1+width, -1, ny, 0, nz-1);
		par->index_ved = i3tensor(start-width, end-1+width, 0, ny, 0, nz-1);
	}
	
	if (par->fd_order == 4){
		par->index_cn = i3tensor(start-1, end+1, -1, ny+1, 0, nz-1);
		par->index_ct = i3tensor(start-width, end-1+width, -2, ny+1, 0, nz-1);
		par->index_ued = i3tensor(start-width, end-1+width, -2, ny+1, 0, nz-1);
		par->index_ved = i3tensor(start-width, end-1+width, -2, ny+2, 0, nz-1);
	}
	
	par->index_les = i3tensor(start-1, end-1+1, 0, ny, 0, nz-1);
	
	par->index_uw_lr = imatrix(-1, ny, 0, nz-1); // freed
	par->index_v_lr = imatrix(0, ny, 0, nz-1); // freed
	
	par->index_tb = imatrix(start, end, 0, nz-1); // freed
	par->kz = dvector(0, nz-1); // freed
	
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				par->index_u[i][j][k] = (i-start)*((ny-1)*nz + (nz-1) + 1) + j*nz + k;
			}
		}
	}
    
	for (i=start-width; i<=end-1+width; ++i){
		for (j=-width_uw; j<=ny-1+width_uw; ++j){
			for (k=0; k<=nz-1; ++k){
				par->index_ued[i][j][k] = (i-start+width)*((ny-1+2*width_uw)*nz + (nz-1) + 1) + (j+width_uw)*nz + k;
			}
		}
	}
    
	
	for (i=start; i<=end-1; ++i){
		for (j=1; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				par->index_v[i][j][k] = (i-start)*((ny-1-1)*nz + (nz-1) + 1) + (j-1)*nz + k;
			}
		}
	}
	
	
	if (par->fd_order == 2){
        for (i=start-width; i<=end-1+width; ++i){
            for (j=1-width_v; j<=ny-1+width_v; ++j){
                for (k=0; k<=nz-1; ++k){
                    par->index_ved[i][j][k] = (i-start+width)*((ny-2+2*width_v)*nz + (nz-1) + 1) + (j-1+width_v)*nz + k;
                    
                }
            }
        }
	}
	
	if (par->fd_order == 4){
        for (i=start-width; i<=end-1+width; ++i){
            for (j=-2; j<=ny+2; ++j){
                for (k=0; k<=nz-1; ++k){
                    par->index_ved[i][j][k] = (i-start+width)*((ny+4)*nz + (nz-1) + 1) + (j+2)*nz + k;
                    
                }
            }
        }
	}
	
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				par->index_wp[i][j][k] = (i-start)*((ny-1)*nz + (nz-1) + 1) + j*nz + k;
			}
		}
	}
	
	for (i=start-width; i<=end-1+width; ++i){
		for (j=-width_uw; j<=ny-1+width_uw; ++j){
			for (k=0; k<=nz-1; ++k){
				par->index_ct[i][j][k] = (i-start+width)*((ny-1+2*width_uw)*nz + (nz-1) + 1) + (j+width_uw)*nz + k;
			}
		}
	}
	
	if (par->fd_order == 2){
		for (i=start; i<=end; ++i){
			for (j=0; j<=ny; ++j){
				for (k=0; k<=nz-1; ++k){
					par->index_cn[i][j][k] = (i-start)*((ny)*nz + (nz-1) + 1) + j*nz + k;
				}
			}
		}
	}
	if (par->fd_order == 4){
		for (i=start-1; i<=end+1; ++i){
			for (j=-1; j<=ny+1; ++j){
				for (k=0; k<=nz-1; ++k){
					par->index_cn[i][j][k] = (i-start+1)*((ny+2)*nz + (nz-1) + 1) + (j+1)*nz + k;
				}
			}
		}
	}
	
	
	for (i=start-1; i<=end; ++i){
		for (j=0; j<=ny; ++j){
			for (k=0; k<=nz-1; ++k){
				par->index_les[i][j][k] = (i-start+1)*((ny)*nz + (nz-1) + 1) + j*nz + k;
			}
		}
	}
	
	for (j=-1; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			par->index_uw_lr[j][k] = (j+1)*(nz) + k;
		}
	}
	for (j=0; j<=ny; ++j){
		for (k=0; k<=nz-1; ++k){
			par->index_v_lr[j][k] = j*(nz) + k;
		}
	}
	
	
	for (i=start; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			par->index_tb[i][k] = (i-start)*(nz) + k;
		}
	}
	
	
	par->kz[0] = 0.0;
	for (k = 1; k<=nz/2; k++ ){
        par->kz[k] = 2*PI*(double)(k)/par->Lz;
		par->kz[nz-k] = 2*PI*(double)(k)/par->Lz;
	}
	
	if (nz%2 == 0){
		par->kz[nz/2] = 0.0;
	}
	
    
	
	par->rk3_alpha[0] =   4.0/15.0;
    par->rk3_alpha[1] =   1.0/15.0;
    par->rk3_alpha[2] =   1.0/ 6.0;
    par->rk3_beta[0]  =   4.0/15.0;
    par->rk3_beta[1]  =   1.0/15.0;
    par->rk3_beta[2]  =   1.0/ 6.0;
	
	
    par->rk3_zeta[0] =        0.0;
    par->rk3_zeta[1] =  -17.0/60.0;
    par->rk3_zeta[2] =   -5.0/12.0;
    par->rk3_gamma[0]  =  8.0/15.0;
    par->rk3_gamma[1]  =  5.0/12.0;
    par->rk3_gamma[2]  =  3.0/ 4.0;
	
	
	if (par->fd_order == 2){
		par->inv_DG = d3tensor(start, end-1, 0, ny-1, 0, nz-1); // freed
	}
	
	if (par->fd_order == 4){
        
		par->inv_DGx0 = d3tensor(start, end-1, 0, ny-1, 0, nz-1);
		par->inv_DGx1 = d3tensor(start, end-1, 0, ny-1, 0, nz-1);
		par->inv_DGx2 = d3tensor(start, end-1, 0, ny-1, 0, nz-1);
        
	}
}


static void readparam(parameters *par)
{
	int i;
	
	if (par->my_rank == 0){
        
		FILE *fp;
		int npar;
        
		fp = fopen("input.par", "r");
		if (fp == NULL){
			printf("parameters_read(): cannot read parameters file");
		}
		
		/************Parameters******************************
         int		nt;				// Total number of time steps
         int		nx;				// # of grids in x-direction
         int		ny;				// # of grids in y-direction
         int		nz;				// # of grids in z-direction
         
         double	Lx;				// size of box in x-direction
         double	Ly;				// size of box in y-direction
         double	Lz;				// size of box in x-direction
         
         double	re;				// Reynolds number
         double	U;				// mean velocity
         double	dt;				// time step size
         double	tol;			// conjugate gradient tolerance
         int		file_dump;		// frecuency of data output
         int		choice;			// output format 1; gnuplot, 2; tecplot, 0; both
         int		const_mass;      // constant mass flux 1; constant pressure gradient 0;
         double	mass_flux;		// consant mass flux
         int		fd_order;		// order of accuracy of the scheme
         ************Parameters*******************************/
        
        
		npar = 0;
		npar += fscanf(fp, "%d    icont",	&par->icont); //add on 16.10.31
		npar += fscanf(fp, "%d    nt",		&par->nt);
		npar += fscanf(fp, "%d   order_of_accuracy",	&par->fd_order);
		npar += fscanf(fp, "%d    nx",		&par->nx);
		npar += fscanf(fp, "%d    ny",		&par->ny);
		npar += fscanf(fp, "%d    nz",		&par->nz);
		npar += fscanf(fp, "%lf   Lx",		&par->Lx);
		npar += fscanf(fp, "%lf   Ly",		&par->Ly);
		npar += fscanf(fp, "%lf   Lz",		&par->Lz);
		npar += fscanf(fp, "%lf   re_aim",		&par->re_aim);
		npar += fscanf(fp, "%lf   re_start",	&par->re_start);
		npar += fscanf(fp, "%lf   U",	    &par->U);
		npar += fscanf(fp, "%lf   dt",	    &par->dt);
		npar += fscanf(fp, "%d    file_dump",	&par->file_dump);
		npar += fscanf(fp, "%d    choice",	&par->choice);
		npar += fscanf(fp, "%d    stat",	&par->stat);
		npar += fscanf(fp, "%lf    dPdx",	&par->dPdx);

		npar += fscanf(fp, "%d    const_mass",	&par->const_mass);

		npar += fscanf(fp, "%lf   mass_flux",	&par->mass_flux_aim);

		npar += fscanf(fp, "%d   les",	&par->les);

		npar += fscanf(fp, "%d   bc",	&par->bc);

		npar += fscanf(fp, "%lf   h0/dy",	&par->h0);

		npar += fscanf(fp, "%d   ode_choice",	&par->ode_choice);

		npar += fscanf(fp, "%d	  j_log",	&par->j_log);

		npar += fscanf(fp, "%d	  bc_data",	&par->bc_data);

		npar += fscanf( fp, "%127s outputdir", par->outputdir );

		npar += fscanf( fp, "%lf epsilon", &par->eps );
        
		if (npar != 27){
			printf("*********************************************\n");
			printf("parameters_read(): failure reading parameters,npar=%d\n",npar);
			printf("*********************************************\n");
		}
		fclose(fp);
        
		/************Parameters******************************
         int		local_size_u;	// order of velocity vector in x-direction
         int		local_size_v;	// order of velocity vector in y-direction
         int		local_size_wp;	// order of velocity vector in w-direction
         ************Parameters*******************************/
		//**********parameters input*************************
		par->local_nx = par->nx/par->nproc;
		par->nq_x = (par->nx-1)*(par->ny)*(par->nz);
		par->nq_y = (par->nx)*(par->ny-1)*(par->nz);
		par->nq_z = (par->nx)*(par->ny)*(par->nz);
		par->np = (par->nx)*(par->ny)*(par->nz);
		par->re = par->U*par->re;
		par->dx = par->Lx/(double)par->nx;
		par->dy = par->Ly/(double)par->ny;
		par->dz = par->Lz/(double)par->nz;
		par->h0 = par->h0*par->dy; //此地已经乘以网格宽度了，正确
		par->re = par->re_start;
		//**********parameters input*************************
        
		printf("\n");
		printf("\n");
		printf("\n");
		printf("icont		        : %d\n",       par->icont);
		printf("nt			: %d\n",       par->nt);
		printf("order of accuracy	: %d\n",       par->fd_order);
		printf("nx			: %d\n",       par->nx);
		printf("ny			: %d\n",       par->ny);
		printf("nz			: %d\n",       par->nz);
		printf("Lx			: %lf\n",      par->Lx);
		printf("Ly			: %lf\n",      par->Ly);
		printf("Lz			: %lf\n",      par->Lz);
		printf("dx			: %lf\n",      par->dx);
		printf("dy			: %lf\n",      par->dy);
		printf("dz			: %lf\n",      par->dz);
		printf("U			: %lf\n",      par->U);
		printf("re_aim			: %lf\n",      par->re_aim);
		printf("re_start		: %lf\n",      par->re_start);
		printf("dt			: %e\n",       par->dt);
		printf("file_dump		: %d\n",       par->file_dump);
		printf("choice			: %d\n",       par->choice);
		printf("stat			: %d\n",       par->stat);
		printf("dPdx			: %lf\n",       par->dPdx);
		printf("mass_flux		: %lf\n",       par->mass_flux_aim);
		printf("const_mass		: %d\n",       par->const_mass);
		printf("les			: %d\n",       par->les);
		printf("bc			: %d\n",       par->bc);
		printf("h0			: %lf\n",       par->h0);
		printf("ode_choice		: %d\n",       par->ode_choice);
		printf("j_log			: %d\n",       par->j_log);
		printf("bc_data			: %d\n",       par->bc_data);
		printf( "outputdir		: %s\n",       par->outputdir );
		printf( "epsilon		: %lf\n",       par->eps );
		
		printf("\n");
		printf("\n");
		printf("\n");
		fflush(stdout);
		
	}
	
    MPI_Bcast( &par->icont,      1, MPI_INT,    0, MPI_COMM_WORLD );	
    MPI_Bcast( &par->nt,         1, MPI_INT,    0, MPI_COMM_WORLD );
	MPI_Bcast( &par->nx,         1, MPI_INT,    0, MPI_COMM_WORLD );
    MPI_Bcast( &par->ny,         1, MPI_INT,    0, MPI_COMM_WORLD );
    MPI_Bcast( &par->nz,         1, MPI_INT,    0, MPI_COMM_WORLD );
    MPI_Bcast( &par->Lx,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &par->Ly,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->Lz,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &par->re_start,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->re_aim,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->re,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->U,		     1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->dt,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->file_dump,  1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->choice,     1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->stat,		 1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->dPdx,		 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->const_mass, 1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->mass_flux_aim,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->fd_order,   1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->les,		 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->bc,		 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->h0,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->ode_choice, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->j_log,      1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->bc_data,    1, MPI_INT, 0, MPI_COMM_WORLD );
	
	MPI_Bcast( &par->local_nx, 1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->nq_x,	   1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->nq_y,	   1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->nq_z,	   1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->np,	   1, MPI_INT   , 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->dx,	   1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->dy,	   1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->dz,	   1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->eps,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
	
	par->mean_kappa = 0.0;
	
	par->width_share = par->fd_order - 1;
	
	int width_v;
	int width_uw;
	const int width = par->width_share;
	
	if (par->fd_order == 2){
		width_v = 1;
		width_uw = 1;
	}
	if (par->fd_order == 4){
		width_v = 3;
		width_uw = 2;
	}
	
	
	par->local_nx_start = par->local_nx*par->my_rank;
	
	par->local_size_u = par->local_nx*par->ny*par->nz;
	par->local_size_v = par->local_nx*(par->ny-1)*par->nz;
	par->local_size_wp = par->local_nx*par->ny*par->nz;
	
	par->local_size_les = (par->local_nx+2)*(par->ny+1)*par->nz;
	
	par->local_size_tb = (par->local_nx+1)*par->nz;
    
	
	// [start-width->end-1+width][-width_uw->ny-1+width_uw][0->nz-1]
	par->local_size_ued = (par->local_nx+2*width)*(par->ny+2*width_uw)*par->nz;
	//par->local_size_ued = (par->local_nx+2)*(par->ny+2)*par->nz;
	
	// [start-width->end-1+width][1-width_v->ny-1+width_v][0->nz-1]
	if (par->fd_order == 2) par->local_size_ved = (par->local_nx+2*width)*(par->ny-1+2*width_v)*par->nz;
	if (par->fd_order == 4) par->local_size_ved = (par->local_nx+2*width)*(par->ny+5)*par->nz;
	
	// [start-width->end-1+width][-width_uw->ny-1+width_uw][0->nz-1]
	par->local_size_ct = (par->local_nx+2*width)*(par->ny+2*width_uw)*par->nz;
	
	// [start->end][0->ny][0->nz-1]
	if (par->fd_order == 2) par->local_size_cn = (par->local_nx+1)*(par->ny+1)*par->nz;
	if (par->fd_order == 4) par->local_size_cn = (par->local_nx+3)*(par->ny+3)*par->nz;
    
	if ( par->my_rank == 0 ) {
        printf( "Process local_nx_start	local_nx	" );
        printf( "local_size_u	local_size_v	local_size_wp\n" );
        fflush( stdout );
    }
    
    for ( i = 0; i <= par->nproc - 1 ; ++i) {
        if ( i == par->my_rank ) {
            printf( "%3d/%3d	%8d	%14d	%8d	%14d	%10d\n",
                   par->my_rank,
                   par->nproc,
                   par->local_nx_start,
                   par->local_nx,
                   par->local_size_u,
                   par->local_size_v,
                   par->local_size_wp );
            fflush( stdout );
        }
        MPI_Barrier( MPI_COMM_WORLD );
    }
    
    if ( par->my_rank == 0 ) {
        printf( "\n" );
        fflush( stdout );
    }
	
	
    if ( par->my_rank == 0 ) {
        par->recvcounts_u = ivector( 0, par->nproc - 1 );
        par->displs_u     = ivector( 0, par->nproc - 1 );
		par->recvcounts_v = ivector( 0, par->nproc - 1 );
        par->displs_v     = ivector( 0, par->nproc - 1 );
		par->recvcounts_wp= ivector( 0, par->nproc - 1 );
        par->displs_wp    = ivector( 0, par->nproc - 1 );
		par->recvcounts_tb= ivector( 0, par->nproc - 1 );
        par->displs_tb    = ivector( 0, par->nproc - 1 );
		
        par->recvcounts_mean_2d= ivector( 0, par->nproc - 1 );
        par->displs_mean_2d    = ivector( 0, par->nproc - 1 );
		par->recvcounts_mean_tb= ivector( 0, par->nproc - 1 );
        par->displs_mean_tb    = ivector( 0, par->nproc - 1 );
    }
    
    MPI_Gather( &par->local_size_u, 1, MPI_INT,
               par->recvcounts_u,   1, MPI_INT,
               0, MPI_COMM_WORLD );
	MPI_Gather( &par->local_size_v, 1, MPI_INT,
               par->recvcounts_v,   1, MPI_INT,
               0, MPI_COMM_WORLD );
	MPI_Gather( &par->local_size_wp, 1, MPI_INT,
               par->recvcounts_wp,   1, MPI_INT,
               0, MPI_COMM_WORLD );
	MPI_Gather( &par->local_size_tb, 1, MPI_INT,
			   par->recvcounts_tb,   1, MPI_INT,
			   0, MPI_COMM_WORLD );
	
	int size;
	size = par->local_nx;
	MPI_Gather( &size, 1, MPI_INT,
			   par->recvcounts_mean_tb,1, MPI_INT,
			   0, MPI_COMM_WORLD );
    size = par->ny*par->local_nx;
	MPI_Gather(&size, 1, MPI_INT,
			   par->recvcounts_mean_2d, 1, MPI_INT,
			   0, MPI_COMM_WORLD );
    
    
    if ( par->my_rank == 0 ) {
        par->displs_u[0] = 0;
		par->displs_v[0] = 0;
		par->displs_wp[0] = 0;
		par->displs_tb[0] = 0;
		
		par->displs_mean_tb[0] = 0;
        par->displs_mean_2d[0] = 0;
        for (i=1; i<=par->nproc-1; ++i){
            par->displs_u[i]  = par->displs_u[i-1]  + par->recvcounts_u[i-1];
			par->displs_v[i]  = par->displs_v[i-1]  + par->recvcounts_v[i-1];
			par->displs_wp[i] = par->displs_wp[i-1] + par->recvcounts_wp[i-1];
			par->displs_tb[i] = par->displs_tb[i-1] + par->recvcounts_tb[i-1];
			
			par->displs_mean_tb[i] = par->displs_mean_tb[i-1] + par->recvcounts_mean_tb[i-1];
            par->displs_mean_2d[i] = par->displs_mean_2d[i-1] + par->recvcounts_mean_2d[i-1];
		}
	}
    
	return;
}



extern void init_vel_bc(variables *var, parameters *par)
{
	int i;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	
	if ( par->my_rank == par->nproc - 1){
		
		var->bc_u_r = FFTW_MALLOC((ny+2)*nz); // freed
		var->bc_v_r = FFTW_MALLOC((ny+1)*(nz)); // freed
		var->bc_w_r = FFTW_MALLOC((ny+2)*nz); // freed
		
		for (i=0; i<=(ny+2)*nz-1; ++i){
			var->bc_u_r[i] = 0.0;
			var->bc_w_r[i] = 0.0;
		}
		
		for (i=0; i<=(ny+1)*nz-1; ++i){
			var->bc_v_r[i] = 0.0;
		}
	}
	
	
	
	if ( par->my_rank == 0){
		var->bc_u_l = FFTW_MALLOC((ny+2)*nz); // freed
		var->bc_v_l = FFTW_MALLOC((ny+1)*(nz)); // freed
		var->bc_w_l = FFTW_MALLOC((ny+2)*nz); // freed
        
		for (i=0; i<=(ny+2)*nz-1; ++i){
			var->bc_u_l[i] = 0.0;
			var->bc_w_l[i] = 0.0;
		}
		
		for (i=0; i<=(ny+1)*nz-1; ++i){
			var->bc_v_l[i] = 0.0;
		}
	}
	
	var->bc_u_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_v_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_w_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_u_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_v_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_w_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	
	var->bc_kappa_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_kappa_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_K_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_K_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_txy_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_txy_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	
	var->bc_rhs_dudy_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_rhs_dudy_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	
	var->bc_dudy_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_dudy_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	
	var->bc_convective_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->bc_convective_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	
	var->Hbc_dudy_t = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	var->Hbc_dudy_b = FFTW_MALLOC((local_nx+1)*(nz)); // freed
	
	
	var->eta0_bar_t = FFTW_MALLOC((local_nx+1)*(nz));
	var->eta0_t = FFTW_MALLOC((local_nx+1)*(nz));
	var->lambda_t = FFTW_MALLOC((local_nx+1)*(nz));
	var->eta0_bar_b = FFTW_MALLOC((local_nx+1)*(nz));
	var->eta0_b = FFTW_MALLOC((local_nx+1)*(nz));
	var->lambda_b = FFTW_MALLOC((local_nx+1)*(nz));
	
	// top and bottom boundary conditions
	for (i=0; i<=(local_nx+1)*(nz)-1; ++i){
		var->bc_u_t[i] = 0.0;
		var->bc_v_t[i] = 0.0;
		var->bc_w_t[i] = 0.0;
		
		var->bc_u_b[i] = 0.0;
		var->bc_v_b[i] = 0.0;
		var->bc_w_b[i] = 0.0;
		
		var->bc_kappa_t[i] = 0.0;
		var->bc_kappa_b[i] = 0.0;
		var->bc_K_t[i] = 0.0;
		var->bc_K_b[i] = 0.0;
		var->bc_txy_t[i] = 0.0;
		var->bc_txy_b[i] = 0.0;
		
		var->bc_dudy_t[i] = 0.0;
		var->bc_dudy_b[i] = 0.0;
		
		var->bc_convective_t[i] = 0.0;
		var->bc_convective_b[i] = 0.0;
		
		var->Hbc_dudy_t[i] = 0.0;
		var->Hbc_dudy_b[i] = 0.0;
		
		var->bc_rhs_dudy_t[i] = 0.0;
		var->bc_rhs_dudy_b[i] = 0.0;
		
		var->eta0_bar_t[i] = 0.0;
		var->eta0_t[i] = 0.0;
		var->lambda_t[i] = 0.0;
		
		var->eta0_bar_b[i] = 0.0;
		var->eta0_b[i] = 0.0;
		var->lambda_b[i] = 0.0;
		
	}
	
	return;
}





extern void init_var(
                     int flag,
                     variables *var,
                     parameters *par)
{
    
	int i;
	const int local_nx = par->local_nx;
	const int local_size_u = par->local_size_u;
	const int local_size_v = par->local_size_v;
	const int local_size_wp = par->local_size_wp;
	const int local_size_cn = par->local_size_cn;
	
	if (flag != 9){
		
		var->u = FFTW_MALLOC(local_size_u); // freed
		var->v = FFTW_MALLOC(local_size_v); // freed
		var->w = FFTW_MALLOC(local_size_wp); // freed
	}
	
	if (flag == 0){
		var->u_phys = FFTW_MALLOC(local_size_u); // freed
		var->v_phys = FFTW_MALLOC(local_size_v); // freed
		var->w_phys = FFTW_MALLOC(local_size_wp); // freed
	}
	
	if (flag != 1){
		var->Nu = FFTW_MALLOC(local_size_u); // freed
		var->Nv = FFTW_MALLOC(local_size_v); // freed
		var->Nw = FFTW_MALLOC(local_size_wp); // freed
	}
    
	
	if (flag == 0){
		var->Lu = FFTW_MALLOC(local_size_u); // freed
		var->Lv = FFTW_MALLOC(local_size_v); // freed
		var->Lw = FFTW_MALLOC(local_size_wp); // freed
		var->bc_Lu = FFTW_MALLOC(local_size_u); // freed
		var->bc_Lv = FFTW_MALLOC(local_size_v); // freed
		var->bc_Lw = FFTW_MALLOC(local_size_wp); // freed
		
		var->Txx = FFTW_MALLOC(local_size_cn);
		var->Tyy = FFTW_MALLOC(local_size_cn);
		var->Tzz = FFTW_MALLOC(local_size_cn);
		var->Txy = FFTW_MALLOC(local_size_cn);
		var->Tyz = FFTW_MALLOC(local_size_cn);
		var->Tzx = FFTW_MALLOC(local_size_cn);
		var->K = FFTW_MALLOC(local_size_cn);
		
		if (par->les == 2){
			var->dudx = FFTW_MALLOC(local_size_wp);
			var->dvdx = FFTW_MALLOC(local_size_wp);
			var->dwdx = FFTW_MALLOC(local_size_wp);
			var->dudy = FFTW_MALLOC(local_size_wp);
			var->dvdy = FFTW_MALLOC(local_size_wp);
			var->dwdy = FFTW_MALLOC(local_size_wp);
			var->dudz = FFTW_MALLOC(local_size_wp);
			var->dvdz = FFTW_MALLOC(local_size_wp);
			var->dwdz = FFTW_MALLOC(local_size_wp);
			
		}
		
		
		if (par->les == 1){
			var->dudx = FFTW_MALLOC(local_size_cn);
			var->dvdx = FFTW_MALLOC(local_size_cn);
			var->dwdx = FFTW_MALLOC(local_size_cn);
			var->dudy = FFTW_MALLOC(local_size_cn);
			var->dvdy = FFTW_MALLOC(local_size_cn);
			var->dwdy = FFTW_MALLOC(local_size_cn);
			var->dudz = FFTW_MALLOC(local_size_cn);
			var->dvdz = FFTW_MALLOC(local_size_cn);
			var->dwdz = FFTW_MALLOC(local_size_cn);
			
		}
		
		var->A1_mean = FFTW_MALLOC(local_nx);
		var->A2_mean = FFTW_MALLOC(local_nx);
		var->B1_mean = FFTW_MALLOC(local_nx);
		var->B2_mean = FFTW_MALLOC(local_nx);
		var->C1_mean = FFTW_MALLOC(local_nx);
		var->C2_mean = FFTW_MALLOC(local_nx);
		var->D_mean = FFTW_MALLOC(local_nx);
		var->E1_mean = FFTW_MALLOC(local_nx);
		var->E2_mean = FFTW_MALLOC(local_nx);
		
	}
	
    
	
	if (flag != 9){
		var->Lu_y = FFTW_MALLOC(local_size_u); // freed
		var->Lv_y = FFTW_MALLOC(local_size_v); // freed
		var->Lw_y = FFTW_MALLOC(local_size_wp); // freed
		var->bc_Lu_y = FFTW_MALLOC(local_size_u); // freed
		var->bc_Lv_y = FFTW_MALLOC(local_size_v); // freed
		var->bc_Lw_y = FFTW_MALLOC(local_size_wp); // freed
	}
	
	if (flag != 1){
		var->Lu_xz = FFTW_MALLOC(local_size_u); // freed
		var->Lv_xz = FFTW_MALLOC(local_size_v); // freed
		var->Lw_xz = FFTW_MALLOC(local_size_wp); // freed
		var->bc_Lu_xz = FFTW_MALLOC(local_size_u); // freed
		var->bc_Lv_xz = FFTW_MALLOC(local_size_v); // freed
		var->bc_Lw_xz = FFTW_MALLOC(local_size_wp); // freed
	}
	
	
	if (flag != 9){
		var->bc_Du = FFTW_MALLOC(local_size_wp); // freed
		var->Du = FFTW_MALLOC(local_size_wp); // freed
        
	}
	
	if (flag != 9){
		var->p = FFTW_MALLOC(local_size_wp); // freed
		var->Gpx = FFTW_MALLOC(local_size_u); // freed
		var->Gpy = FFTW_MALLOC(local_size_v); // freed
		var->Gpz = FFTW_MALLOC(local_size_wp); // freed
		var->DGp = FFTW_MALLOC(local_size_wp); // freed
	}
	
	
	if (flag == 0){
		for (i=0; i<=local_size_u-1; ++i){
			var->u[i] = 0.0;
			var->u_phys[i] = 0.0;
			var->Nu[i] = 0.0;
			var->Lu[i] = 0.0;
			var->bc_Lu[i] = 0.0;
			var->Lu_y[i] = 0.0;
			var->bc_Lu_y[i] = 0.0;
			var->Lu_xz[i] = 0.0;
			var->bc_Lu_xz[i] = 0.0;
			var->Gpx[i] = 0.0;
		}
        
		for (i=0; i<=local_size_v-1; ++i){
			var->v[i] = 0.0;
			var->v_phys[i] = 0.0;
			var->Nv[i] = 0.0;
			var->Lv[i] = 0.0;
			var->bc_Lv[i] = 0.0;
			var->Lv_y[i] = 0.0;
			var->bc_Lv_y[i] = 0.0;
			var->Lv_xz[i] = 0.0;
			var->bc_Lv_xz[i] = 0.0;
			var->Gpy[i] = 0.0;
		}
        
		for (i=0; i<=local_size_wp-1; ++i){
			var->w[i] = 0.0;
			var->w_phys[i] = 0.0;
			var->Nw[i] = 0.0;
			var->Lw[i] = 0.0;
			var->bc_Lw[i] = 0.0;
			var->Lw_y[i] = 0.0;
			var->bc_Lw_y[i] = 0.0;
			var->Lw_xz[i] = 0.0;
			var->bc_Lw_xz[i] = 0.0;
            
			var->bc_Du[i] = 0.0;
			var->Du[i] = 0.0;
            
			var->p[i] = 0.0;
            
			var->DGp[i] = 0.0;
			var->Gpz[i] = 0.0;
		}
		
		for (i=0; i<=local_size_cn-1; ++i){
			var->Txx[i] = 0.0;
			var->Tyy[i] = 0.0;
			var->Tzz[i] = 0.0;
			
			var->Txy[i] = 0.0;
			var->Tyz[i] = 0.0;
			var->Tzx[i] = 0.0;
			
			var->K[i] = 0.0;
		}
		
		
		if (par->les == 2){
			for (i=0; i<=local_size_wp-1; ++i){
				var->dudx[i] = 0.0;
				var->dvdx[i] = 0.0;
				var->dwdx[i] = 0.0;
				var->dudy[i] = 0.0;
				var->dvdy[i] = 0.0;
				var->dwdy[i] = 0.0;
				var->dudz[i] = 0.0;
				var->dvdz[i] = 0.0;
				var->dwdz[i] = 0.0;
			}
		}
		
		if (par->les == 1){
			for (i=0; i<=local_size_cn-1; ++i){
				var->dudx[i] = 0.0;
				var->dvdx[i] = 0.0;
				var->dwdx[i] = 0.0;
				var->dudy[i] = 0.0;
				var->dvdy[i] = 0.0;
				var->dwdy[i] = 0.0;
				var->dudz[i] = 0.0;
				var->dvdz[i] = 0.0;
				var->dwdz[i] = 0.0;
			}
		}
		
		
		
	}
	
	
	if (flag == 1){
		for (i=0; i<=local_size_u-1; ++i){
			var->u[i] = 0.0;
			var->Lu_y[i] = 0.0;
			var->bc_Lu_y[i] = 0.0;
			var->Gpx[i] = 0.0;
		}
        
		for (i=0; i<=local_size_v-1; ++i){
			var->v[i] = 0.0;
			var->Lv_y[i] = 0.0;
			var->bc_Lv_y[i] = 0.0;
			var->Gpy[i] = 0.0;
		}
        
		for (i=0; i<=local_size_wp-1; ++i){
			var->w[i] = 0.0;
			var->Lw_y[i] = 0.0;
			var->bc_Lw_y[i] = 0.0;
            
			var->bc_Du[i] = 0.0;
			var->Du[i] = 0.0;
            
			var->p[i] = 0.0;
            
			var->DGp[i] = 0.0;
			var->Gpz[i] = 0.0;
		}
	}
	
	
	if (flag == 9){
		for (i=0; i<=local_size_u-1; ++i){
			var->Nu[i] = 0.0;
			var->Lu_xz[i] = 0.0;
			var->bc_Lu_xz[i] = 0.0;
		}
        
		for (i=0; i<=local_size_v-1; ++i){
			var->Nv[i] = 0.0;
			var->Lv_xz[i] = 0.0;
			var->bc_Lv_xz[i] = 0.0;
		}
        
		for (i=0; i<=local_size_wp-1; ++i){
			var->Nw[i] = 0.0;
			var->Lw_xz[i] = 0.0;
			var->bc_Lw_xz[i] = 0.0;
		}
	}
	
    
}



extern void finalize_var(
                         int flag,
                         variables *var,
                         parameters *par)
{
	
	if (flag !=9){
		fftw_free(var->u);
		fftw_free(var->v);
		fftw_free(var->w);
	}
    
	if (flag ==0){
		fftw_free(var->u_phys);
		fftw_free(var->v_phys);
		fftw_free(var->w_phys);
	}
	
	if (flag !=1){
		fftw_free(var->Nu);
		fftw_free(var->Nv);
		fftw_free(var->Nw);
	}
	
	if (flag == 0){
		fftw_free(var->Lu);
		fftw_free(var->Lv);
		fftw_free(var->Lw);
		fftw_free(var->bc_Lu);
		fftw_free(var->bc_Lv);
		fftw_free(var->bc_Lw);
		
		fftw_free(var->Txx);
		fftw_free(var->Tyy);
		fftw_free(var->Tzz);
		fftw_free(var->Txy);
		fftw_free(var->Tyz);
		fftw_free(var->Tzx);
		
		fftw_free(var->K);
		
		fftw_free(var->A1_mean);
		fftw_free(var->A2_mean);
		fftw_free(var->B1_mean);
		fftw_free(var->B2_mean);
		fftw_free(var->C1_mean);
		fftw_free(var->C2_mean);
		fftw_free(var->D_mean);
		fftw_free(var->E1_mean);
		fftw_free(var->E2_mean);
		
		if (par->les != 0){
			fftw_free(var->dudx);
			fftw_free(var->dvdx);
			fftw_free(var->dwdx);
			fftw_free(var->dudy);
			fftw_free(var->dvdy);
			fftw_free(var->dwdy);
			fftw_free(var->dudz);
			fftw_free(var->dvdz);
			fftw_free(var->dwdz);
		}
        
	}
	
	if (flag !=9){
		fftw_free(var->Lu_y);
		fftw_free(var->Lv_y);
		fftw_free(var->Lw_y);
		fftw_free(var->bc_Lu_y);
		fftw_free(var->bc_Lv_y);
		fftw_free(var->bc_Lw_y);
	}
	
	if (flag !=1){
		fftw_free(var->Lu_xz);
		fftw_free(var->Lv_xz);
		fftw_free(var->Lw_xz);
		fftw_free(var->bc_Lu_xz);
		fftw_free(var->bc_Lv_xz);
		fftw_free(var->bc_Lw_xz);
	}
	
	if (flag !=9){
		fftw_free(var->bc_Du);
		fftw_free(var->Du);
	}
	
	if (flag !=9){
		fftw_free(var->p);
		fftw_free(var->Gpx);
		fftw_free(var->Gpy);
		fftw_free(var->Gpz);
		fftw_free(var->DGp);
	}
	
}


extern void finalize_par(
                         parameters *par)
{
	const int local_nx_start = par->local_nx_start;
	const int local_nx = par->local_nx;
	const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	const int width = par->width_share;
	int width_v;
	int width_uw;
	
	if (par->fd_order == 2){
		width_v = 1;
		width_uw = 1;
	}
	if (par->fd_order == 4){
		width_v = 3;
		width_uw = 2;
	}
	
	free_i3tensor(par->index_u, start, end-1, 0, ny-1, 0, nz-1);
	free_i3tensor(par->index_v, start, end-1, 1, ny-1, 0, nz-1);
	free_i3tensor(par->index_wp, start, end-1, 0, ny-1, 0, nz-1);
	
	free_i3tensor(par->index_les, start-1, end, 0, ny, 0, nz-1);
	
	free_imatrix(par->index_uw_lr, -1, ny, 0, nz-1);
	free_imatrix(par->index_v_lr, 0, ny, 0, nz-1);
	
	free_imatrix(par->index_tb, start, end, 0, nz-1);
	free_dvector(par->kz, 0, nz-1);
	
	free_i3tensor(par->index_ct,  start-width, end-1+width, -width_uw, ny-1+width_uw, 0, nz-1);
	free_i3tensor(par->index_ued, start-width, end-1+width, -width_uw, ny-1+width_uw, 0, nz-1);
	
	if (par->fd_order == 2){
		free_i3tensor(par->index_ved, start-width, end-1+width, 1-width_v, ny-1+width_v, 0, nz-1);
		free_i3tensor(par->index_cn, start, end, 0, ny, 0, nz-1);
		
		free_d3tensor(par->inv_DG, start, end-1, 0, ny-1, 0, nz-1);
		
	}
	if (par->fd_order == 4){
		free_i3tensor(par->index_ved, start-width, end-1+width, -2, ny+2, 0, nz-1);
		free_i3tensor(par->index_cn, start-1, end+1, -1, ny+1, 0, nz-1);
		
		free_d3tensor(par->inv_DGx0, start, end-1, 0, ny-1, 0, nz-1);
		free_d3tensor(par->inv_DGx1, start, end-1, 0, ny-1, 0, nz-1);
		free_d3tensor(par->inv_DGx2, start, end-1, 0, ny-1, 0, nz-1);
	}
    
	
	if ( par->my_rank == 0 ) {
        free_ivector( par->recvcounts_u, 0, par->nproc - 1 );
        free_ivector( par->displs_u,     0, par->nproc - 1 );
		free_ivector( par->recvcounts_v, 0, par->nproc - 1 );
        free_ivector( par->displs_v,     0, par->nproc - 1 );
		free_ivector( par->recvcounts_wp,0, par->nproc - 1 );
        free_ivector( par->displs_wp,    0, par->nproc - 1 );
        
        free_ivector( par->recvcounts_mean_tb,0, par->nproc - 1 );
        free_ivector( par->displs_mean_tb,    0, par->nproc - 1 );
        free_ivector( par->recvcounts_mean_2d,0, par->nproc - 1 );
        free_ivector( par->displs_mean_2d,    0, par->nproc - 1 );
        
		free_i3tensor(par->index_bigwp, 0, nx-1, 0, ny-1, 0, nz-1);
	}
	
	return;
}


void finalize_vel_bc(
                     variables *var,
                     parameters *par)
{
	
	if ( par->my_rank == 0){
		fftw_free(var->bc_u_l);
		fftw_free(var->bc_v_l);
		fftw_free(var->bc_w_l);
	}
	
	if ( par->my_rank == par->nproc - 1){
		fftw_free(var->bc_u_r);
		fftw_free(var->bc_v_r);
		fftw_free(var->bc_w_r);
	}
    
	fftw_free(var->bc_u_t);
	fftw_free(var->bc_v_t);
	fftw_free(var->bc_w_t);
	fftw_free(var->bc_u_b);
	fftw_free(var->bc_v_b);
	fftw_free(var->bc_w_b);
	
	fftw_free(var->bc_kappa_t);
	fftw_free(var->bc_kappa_b);
	fftw_free(var->bc_K_t);
	fftw_free(var->bc_K_b);
	fftw_free(var->bc_txy_t);
	fftw_free(var->bc_txy_b);
	
	fftw_free(var->bc_dudy_t);
	fftw_free(var->bc_dudy_b);
	
	fftw_free(var->bc_convective_t);
	fftw_free(var->bc_convective_b);
	
	fftw_free(var->Hbc_dudy_t);
	fftw_free(var->Hbc_dudy_b);
	
	fftw_free(var->bc_rhs_dudy_t);
	fftw_free(var->bc_rhs_dudy_b);
	
	fftw_free(var->eta0_bar_t);
	fftw_free(var->eta0_t);
	fftw_free(var->lambda_t);
	
	fftw_free(var->eta0_bar_b);
	fftw_free(var->eta0_b);
	fftw_free(var->lambda_b);
	
	return;
    
}


