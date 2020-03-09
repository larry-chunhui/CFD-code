/*
 * main.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * main() function comprises three major parts.
 * (1) Simulation initializations -- reading constants, allocating necessary memory for the field variables and initializing the flow field variables and parameters.
 *                                -- MPI set-up is also done here.
 * (2) Time stepping -- Please refer to the accompanying document for the overview of time-stepping methods.
 *                   -- Outputting some data. The functions for outputting the data are also defined in this file.
 * (3) Finalizing the simulations -- deallocating the memory.
 * 
 * Time for each iteration is recorded in time_per_iteration.dat. In addition, the total time is printed at the end of the simulation.
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "share.h"
#include "interpolation.h"
#include "interpolation__2.h"
#include "interpolation__4.h"
#include "derivatives__2.h"
#include "derivatives__4.h"
#include "time_march.h"
#include "error_check.h"
#include "profiles.h"
#include "bc__get.h"
#include "velgrad_tensor.h"
#include "les.h"
#include "les_diagnostics.h"
#include "statistics.h"
#include "profiles.h"
#include "spiral.h"

static void time_per_iteration(
                               clock_t *t_now,
                               parameters *par);

void initialize_flow_field(
						   variables *var0,
						   variables *var1,
						   parameters *par,
						   fftwplans *planptr);

void initialize_bc_flow_field(
							  variables *var,
							  parameters *par,
							  fftwplans *planptr);

void output(variables *var, parameters *par);

extern void dump_physical_vector(
                                 double *u_phys,
                                 double *v_phys,
                                 double *w_phys,
                                 parameters *par );


extern void read_physical_vector(      // read initial velocity field
                                 double *u_phys,
                                 double *v_phys,
                                 double *w_phys,
                                 parameters *par );
extern void initialize_velocity_field(  // generate initial random velocity field by code
                                 double *u_phys,
                                 double *v_phys,
                                 double *w_phys,
                                 parameters *par );





extern void bc_dump_physical_vector(
                                    double *bc_u_phys,
                                    double *bc_v_phys,
                                    double *bc_w_phys,
                                    const char *basename,
                                    parameters *par );

extern void bc_u_tau_dump_physical_vector(
                                          double *bc_u_tau,
                                          const char *basename,
                                          parameters *par );


extern void bc_read_physical_vector(
                                    double *bc_u_phys,
                                    double *bc_v_phys,
                                    double *bc_w_phys,
                                    const char *basename,
                                    parameters *par );


void check_divergence(
                      variables *var,
                      parameters *par,
                      fftwplans *planptr);


void channel_integration(
                         variables *var9,
                         variables *var0,
                         variables *var1,
                         parameters *par,
                         fftwplans *planptr);

extern void get_mass_flux(
						  double *u,
						  parameters *par );





extern void dump_plane_data_plot(
                                 double *plane,
                                 const char *basename,
                                 parameters *par);

extern void dump_physical_vector_plot(
									  double *u_phys,
									  double *v_phys,
									  double *w_phys,
									  parameters *par );

/* ******************************************************************* */
/* Main.                                                               */
/* ******************************************************************* */
int main(int argc, char **argv)
{
	time_t t_0;
	clock_t t_now;
    
    
	char hostname[128];
	size_t len = 126;
	pid_t pid;
	
	parameters p, *par;
	variables v9, v0, v1, *var9, *var0, *var1;
	fftwplans	plan;
	fftwplans	*planptr;
	planptr = &plan;
	int	nx, ny, nz;
    
    
	double MPIt1, MPIt2, MPIelapsed;
	
	MPI_Init(&argc, &argv);
	MPIt1 = MPI_Wtime();
	
	t_0 = time(NULL);  // Set the start time.
	t_now = time(NULL);
	
	par = &p;
	var9 = &v9;
	var0 = &v0;
	var1 = &v1;
	
	init_par(par);
	
	init_var(9, var9, par);
	init_var(0, var0, par);
	init_var(1, var1, par);
	
	init_vel_bc(var9, par);
	init_vel_bc(var0, par);
	init_vel_bc(var1, par);
	
	
	fft_make_plans(var0, par, planptr);
	fft_make_plans__transpose(par, planptr);
    
	nx = par->nx;
	ny = par->ny;
	nz = par->nz;
	
	gethostname(hostname, len);
	pid = getpid();
    
    
	/* Print out the hostname and process id for this process */
//	printf("Greeting from process %d(%s:%d) out of %d!\n",par->my_rank,hostname,pid,par->nproc);
	fflush(stdout);
	
	initialize_flow_field(var0, var1, par, planptr);
if(par->my_rank == 0){
//printf("initialize_flow_field succeed\n");
}
	get_mass_flux(var0->u, par);
if(par->my_rank == 0){
//printf("get_mass_flux succeed\n");
}	
	MPI_Barrier(MPI_COMM_WORLD);
    
	par->it = 0;
	par->tott = 0.0;
if(par->my_rank == 0){
//printf("channel_integration begin succeed\n");
}	
    while ( par->it <= par->nt ) {
        
		channel_integration(var9, var0, var1, par, planptr);   //
//if(par->my_rank == 0){
//printf("channel_integration  succeed %d\n",par->it);
//}       
		if (par->my_rank == 0){
			time_per_iteration(&t_now, par);
		}
		par->it += 1;
		par->tott += par->dt;
        
		if (par->it%par->file_dump == 0){
			dump_physical_vector(var0->u_phys, var0->v_phys, var0->w_phys, par);
		}
		if (par->it%par->file_dump == 0 && par->bc == 1){
			dump_physical_vector(var0->u_phys, var0->v_phys, var0->w_phys, par);
			fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_u_t, var9->bc_u_t);
			fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_v_t, var9->bc_v_t);
			fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_w_t, var9->bc_w_t);
			fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_u_b, var9->bc_u_b);
			fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_v_b, var9->bc_v_b);
			fftw_execute_r2r(planptr->p1d_invz_bc_tb, var1->bc_w_b, var9->bc_w_b);
			bc_dump_physical_vector(var9->bc_u_t, var9->bc_v_t, var9->bc_w_t, "top", par);
			bc_dump_physical_vector(var9->bc_u_b, var9->bc_v_b, var9->bc_w_b, "bot", par);
			bc_u_tau_dump_physical_vector(var0->bc_dudy_t, "top", par);
			bc_u_tau_dump_physical_vector(var0->bc_dudy_b, "bot", par);
		}
	}
	
	
	if (par->my_rank == 0){
		output(var0, par);
	}
	
	finalize_vel_bc(var9, par);
	finalize_vel_bc(var0, par);
	finalize_vel_bc(var1, par);
	finalize_var(9, var9, par);
	finalize_var(0, var0, par);
	finalize_var(1, var1, par);
	
	fft_destroy_plans(planptr, par);
    	fft_destroy_plans__transpose(planptr, par);
	finalize_par(par);
	
	
	//print out time passed
	MPIt2 = MPI_Wtime();
	MPIelapsed=MPIt2-MPIt1;
    
	if( par->my_rank == 0 ) {
		printf("\ttime MPI     = %f s\n",MPIelapsed);
	}
	
	MPI_Finalize();
	exit(0);
	
}





void channel_integration(
                         variables *var9,
                         variables *var0,
                         variables *var1,
                         parameters *par,
                         fftwplans *planptr)
{
	
	int i;
	double *u_copy, *v_copy, *w_copy;
	double *res_diss, *sub_diss;
	
	// Control for mass flux
	const double trans_T = 10.0; // Transisiton time
    
	par->re = par->re_aim
    + (par->re_start - par->re_aim)*exp(-4.0*par->tott/trans_T);
	
	
	if (par->tott <= trans_T) {
		par->mass_flux = par->mass_flux_start
        + (par->mass_flux_aim - par->mass_flux_start)/trans_T*par->tott;
    } else {
        par->mass_flux = par->mass_flux_aim;
    }
    
    
	
	
	if (par->my_rank == 0){
		printf("time step at %d dt = %e T = %e\n", \
			   par->it, par->dt, par->tott);
		printf("Re_aim = %e Re_curr = %e\n", \
			   par->re_aim, par->re);
		printf("flux_aim = %e flux_curr = %e\n", \
			   par->mass_flux_aim, par->mass_flux);
		printf("\n");
	}
// 已经成功	
	if (par->bc == 1 && par->my_rank == 0){
		printf("wall ode terms on average\n");
		printf("A	%e	%e\n", par->A1, par->A2);
		printf("B	%e	%e\n", par->B1, par->B2);
		printf("C	%e	%e\n", par->C1, par->C2);
		printf("D	%e\n", par->D);
		printf("E	%e	%e\n", par->E1, par->E2);
		printf("\n");
	}
	
// if(par->my_rank == 0){
//	printf(" channel integration /begin time march 0");
//}	
	time_march(var9, var0, var1, par, 0, planptr);
	MPI_Barrier(MPI_COMM_WORLD);
 //if(par->my_rank == 0){
//	printf(" channel integration /begin time march 1");
//}		
	time_march(var9, var0, var1, par, 1, planptr);
	MPI_Barrier(MPI_COMM_WORLD);

// if(par->my_rank == 0){
//	printf(" channel integration /begin time march 2");
//}		
	time_march(var9, var0, var1, par, 2, planptr);
	MPI_Barrier(MPI_COMM_WORLD);
 
// if(par->my_rank == 0){
//	printf(" channel integration /end time march 2");
//}	   
	check_divergence(var0, par, planptr);
 
// if(par->my_rank == 0){
//	printf(" channel integration /check_divergence");
//}	
   
	u_copy = FFTW_MALLOC(par->local_size_u); // freed
	v_copy = FFTW_MALLOC(par->local_size_v); // freed
	w_copy = FFTW_MALLOC(par->local_size_wp); // freed
    
	for (i=0; i<=par->local_size_u-1; ++i){
		u_copy[i] = var0->u[i];
	}
	for (i=0; i<=par->local_size_v-1; ++i){
		v_copy[i] = var0->v[i];
	}
	for (i=0; i<=par->local_size_wp-1; ++i){
		w_copy[i] = var0->w[i];
	}
    
	fftw_execute_r2r(planptr->p1d_invz_u, u_copy, var0->u_phys);
	fftw_execute_r2r(planptr->p1d_invz_v, v_copy, var0->v_phys);
	fftw_execute_r2r(planptr->p1d_invz_wp, w_copy, var0->w_phys);
    
	fftw_free(u_copy);
	fftw_free(v_copy);
	fftw_free(w_copy);
	MPI_Barrier(MPI_COMM_WORLD);
    
    
	if (par->my_rank == 0 && par->it%par->stat == 0){
		output(var0, par);
	}
    
    	if (par->it%par->stat == 0){
        dump_physical_vector_plot(var0->u_phys, var0->v_phys, var0->w_phys, par);
	}
    
	if (par->it%par->stat == 0 && par->les == 1){
        	res_diss = FFTW_MALLOC(par->local_size_cn); // freed
		sub_diss = FFTW_MALLOC(par->local_size_cn); // freed
		
		fftw_execute_r2r(planptr->p1d_invz_cn, var0->Txx, var0->Txx);
		fftw_execute_r2r(planptr->p1d_invz_cn, var0->Tyy, var0->Tyy);
		fftw_execute_r2r(planptr->p1d_invz_cn, var0->Tzz, var0->Tzz);
		fftw_execute_r2r(planptr->p1d_invz_cn, var0->Txy, var0->Txy);
		fftw_execute_r2r(planptr->p1d_invz_cn, var0->Tyz, var0->Tyz);
		fftw_execute_r2r(planptr->p1d_invz_cn, var0->Tzx, var0->Tzx);
		
		les_get_dissipation(
                            var0->dudx, var0->dudy, var0->dudz,
                            var0->dvdx, var0->dvdy, var0->dvdz,
                            var0->dwdx, var0->dwdy, var0->dwdz,
                            var0->Txx, var0->Tyy, var0->Tzz,
                            var0->Txy, var0->Tyz, var0->Tzx,
                            res_diss, sub_diss, par);
        
		mean_profile_eta(res_diss, par);
		mean_dissipation_ratio(res_diss, sub_diss, par);
		
		fftw_free(res_diss);
		fftw_free(sub_diss);
	}
	
if(par->my_rank == 0){
//printf("skin_friction begin\n");
}    
	
	if (par->bc == 1){
		skin_friction__bc(var0, par);
		dump_bc_statistics(var0, par);
	}else{
		skin_friction(var0->u, par);
		
	}
if(par->my_rank == 0){
// printf("statistics_u begin\n");
}
	statistics_u(var0->u, var0->v, var0->w, par);
	MPI_Barrier(MPI_COMM_WORLD);
    
	
	courant(var0->u_phys, var0->v_phys, var0->w_phys, par);
	
	// Control for cfl.
	const double omega_cfl = 0.7; // Under-relaxation
	const double want_cfl  = 1.0; // Desired Courant number
	par->dt *= 1.0 + omega_cfl*( want_cfl/par->cfl - 1.0 );
	
	
	if (par->my_rank == 0){
		printf("u_tau = %e CFL = %e\n", \
			   par->u_tau, par->cfl);
		printf("expected dpdx_t(b) = %e	%e	dPdx = %e\n", \
			   par->expected_dpdx_t, par->expected_dpdx_b, par->dPdx);
		printf("bc_u_b = %e bc_u_t = %e	kappa = %e\n", \
			   var0->bc_u_b[par->index_tb[0][0]], var0->bc_u_t[par->index_tb[0][0]], par->mean_kappa);
		printf("bc_txy_t(b) = (%e	%e)	bc_K_t(b) = (%e	%e)\n", \
			   par->mean_bctxy_t, par->mean_bctxy_b, par->mean_bcK_t, par->mean_bcK_b);
		printf("dudy = %e	div = %e\n", \
			   par->mean_dudy, par->Du_max);
		printf("\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
    
	return;
}





static void time_per_iteration(
                               clock_t *t_now,
                               parameters *par)
{
	clock_t t_old = *t_now;
    *t_now = clock();
    
	FILE *fp;
	char s[100];
    
	sprintf(s, "./outputdir/time_per_iteration.dat");
	fp = fopen(s, "a" );
    
	fprintf( fp, "%d %+e\n", par->it,
            (double)( *t_now - t_old )/CLOCKS_PER_SEC );
    
	fclose( fp );
}



void output(variables *var, parameters *par)
{
	
	int i, j;
	int ny = par->ny;
	int start = par->local_nx_start;
	int end = par->local_nx_start + par->local_nx;
	
	char s[100];
	FILE *out;
	
    
	sprintf(s, "./outputdir/u_field_%d.txt", par->it);
	out = fopen(s, "w");
    
	for (i=start; i<=end-1; i++){
		for (j=0; j<=ny-1; j++){
			fprintf(out, "%e	%e	%e\n", i*par->dx, (j+0.5)*par->dy, var->u_phys[par->index_u[i][j][0]]);
            
		}
		fprintf(out, "\n");
	}
	fclose(out);
	
	
    
	sprintf(s, "./outputdir/v_field_%d.txt",par->it);
	out = fopen(s, "w");
    
	for (i=start; i<=end-1; i++){
		for (j=1; j<=ny-1; j++){
			fprintf(out, "%e	%e	%e\n", (i+0.5)*par->dx, j*par->dy, var->v_phys[par->index_v[i][j][0]]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	
	
    
	sprintf(s, "./outputdir/w_field_%d.txt",par->it);
	out = fopen(s, "w");
    
	for (i=start; i<=end-1; i++){
		for (j=0; j<=ny-1; j++){
			fprintf(out, "%e	%e	%e\n", (i+0.5)*par->dx, (j+0.5)*par->dy, var->w_phys[par->index_wp[i][j][0]]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
    
	return;
}


extern void dump_physical_vector(
                                 double *u_phys,
                                 double *v_phys,
                                 double *w_phys,
                                 parameters *par )
{
	int total_size_u;
	int local_size_u = par->local_size_u;
	int total_size_v;
	int local_size_v = par->local_size_v;
	int total_size_wp;
	int local_size_wp = par->local_size_wp;
    double *u_big, *v_big, *w_big;
    
    int nx = par->nx;
    int ny = par->ny;
    int nz = par->nz;
    
	total_size_u = nx*ny*nz;
	total_size_v = nx*(ny-1)*nz;
	total_size_wp = nx*ny*nz;
//if(par->my_rank == 0){
//printf("begin dump_physical_vector succeed\n");
//}
    if ( par->my_rank == 0 ) {
        u_big = (double *)fftw_malloc(sizeof(double)*total_size_u);//freed
        v_big = (double *)fftw_malloc(sizeof(double)*total_size_v);//freed
        w_big = (double *)fftw_malloc(sizeof(double)*total_size_wp);//freed
    }
    else {
        u_big = NULL;
        v_big = NULL;
        w_big = NULL;
    }
    
    MPI_Gatherv( u_phys, local_size_u,             MPI_DOUBLE,
                u_big,  par->recvcounts_u, par->displs_u, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
    MPI_Gatherv( v_phys, local_size_v,             MPI_DOUBLE,
                v_big,  par->recvcounts_v, par->displs_v, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
    MPI_Gatherv( w_phys, local_size_wp,             MPI_DOUBLE,
                w_big,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
    
    if ( par->my_rank == 0 ) {
        
        FILE *fp;
        char filename[100];
        int i, j, k, index;


//printf("begin dump_u_vector succeed\n");
//printf("%s\n",par->outputdir);
/*		sprintf( filename, "./input/re%d/u_dump_%d%d%d.dat", (int)par->re_start, nx, ny, nz);
        fp = fopen( filename, "r" );
        if ( fp == NULL ) printf( "read_physical_vector(): cannot read from file" );
		
		
        fscanf( fp, "%d,%d,%d\n", &i, &j, &k );
        if ( !( i == nx && j == ny && k == nz ) ){
           printf("U!! i=%d, j=%d, k=%d, but expecting %d, %d, %d\n", i,j,k,nx,ny , nz);
			printf( "read_physical_vector(): inconsistent dimensions" );
		}
		
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					fscanf( fp, "%lf\n", &u_big[index]);
					index = index + 1;
                }
			}
		}
		
        fclose( fp );*/    
        
//	 sprintf( filename, "./%s/u_dump_it%d.dat", par->outputdir, par->it );
	 sprintf( filename, "./outputdir/u_dump_it%d.dat",par->it );  //I change 9.26

        fp = fopen( filename, "w" );

// printf("fopen dump_physical_vector succeed\n");
      
        fprintf( fp, "%d,%d,%d\n", nx, ny, nz);
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
                    fprintf( fp, "%lf\n",   u_big[index]);
					index = index + 1;
				}
			}
		}
        
        fclose( fp );


//printf("begin dump_v_vector succeed\n");
    		
        
//		sprintf( filename, "./%s/v_dump_it%d.dat", par->outputdir, par->it );
		sprintf( filename, "./outputdir/v_dump_it%d.dat",par->it );  //I change 9.26

        fp = fopen( filename, "w" );
        
        
        fprintf( fp, "%d,%d,%d\n", nx, ny-1, nz);
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
                    fprintf( fp, "%lf\n",   v_big[index]);
					index = index + 1;
				}
			}
		}
        
        fclose( fp );

//if(par->my_rank == 0){
//printf("begin dump_u_vector succeed\n");
//}     		
        
//		sprintf( filename, "./%s/w_dump_it%d.dat", par->outputdir, par->it );
		sprintf( filename, "./outputdir/w_dump_it%d.dat",par->it );  //I change 9.26
        fp = fopen( filename, "w" );
        
        
        fprintf( fp, "%d,%d,%d\n", nx, ny, nz);
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
                    fprintf( fp, "%lf\n",   w_big[index]);
					index = index + 1;
				}
			}
		}
        
        fclose( fp );
    }
    
    if ( par->my_rank == 0 ) {
        fftw_free(u_big);
        fftw_free(v_big);
        fftw_free(w_big);
    }
	
	return;
}






extern void read_physical_vector(
                                 double *u_phys,
                                 double *v_phys,
                                 double *w_phys,
                                 parameters *par )
{
    int nx, ny, nz,  icont;
	int total_size_u;
	int total_size_v;
	int total_size_wp;
    double *u_big, *v_big, *w_big;
    
    nx = par->nx;
    ny = par->ny;
    nz = par->nz;
	icont = par->icont;
    
	total_size_u = nx*ny*nz;
	total_size_v = nx*(ny-1)*nz;
	total_size_wp = nx*ny*nz;
    
	printf( "read_physical_vector(): test\n" );
	
	
    if ( par->my_rank == 0 ) {
        u_big = (double *)fftw_malloc(sizeof(double)*total_size_u);
        v_big = (double *)fftw_malloc(sizeof(double)*total_size_v);
        w_big = (double *)fftw_malloc(sizeof(double)*total_size_wp);
    }
    else {
        u_big = NULL;
        v_big = NULL;
        w_big = NULL;
    }

    if ( par->my_rank == 0 && icont == 1 ) {
        
        FILE *fp;
        char filename[100];
        int i, j, k, index;
		
		sprintf( filename, "./input/re%d/u_dump_%d%d%d.dat", (int)par->re_start, nx, ny, nz);
        fp = fopen( filename, "r" );
        if ( fp == NULL ) printf( "read_physical_vector(): cannot read from file" );
		
		
        fscanf( fp, "%d,%d,%d\n", &i, &j, &k );
        if ( !( i == nx && j == ny && k == nz ) ){
           printf("U!! i=%d, j=%d, k=%d, but expecting %d, %d, %d\n", i,j,k,nx,ny , nz);
			printf( "read_physical_vector(): inconsistent dimensions" );
		}
		
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					fscanf( fp, "%lf\n", &u_big[index]);
					index = index + 1;
                }
			}
		}
		
        fclose( fp );
        
		sprintf( filename, "./input/re%d/v_dump_%d%d%d.dat", (int)par->re_start, nx, ny, nz);
		fp = fopen( filename, "r" );
        	if ( fp == NULL ) printf( "read_physical_vector(): cannot read from file" );
		
		
        fscanf( fp, "%d,%d,%d\n", &i, &j, &k );
		if ( !( i == nx && j == ny-1 && k == nz ) ){
                    printf("i=%d, j=%d, k=%d, but expecting %d, %d, %d\n", i,j,k, nx,ny-1,nz);
			printf( "read_physical_vector(): inconsistent dimensions" );
		}
		
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=1; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					fscanf( fp, "%lf\n", &v_big[index]);
            		index = index + 1;
                }
			}
		}
        
        fclose( fp );
		
	 sprintf( filename, "./input/re%d/w_dump_%d%d%d.dat", (int)par->re_start, nx, ny, nz);
        fp = fopen( filename, "r" );
        if ( fp == NULL ) printf( "read_physical_vector(): cannot read from file" );
		
		
		fscanf( fp, "%d,%d,%d\n", &i, &j, &k );
		if ( !( i == nx && j == ny && k == nz ) ){
            
            		printf("W!!! i=%d, j=%d, k=%d, but expecting %d, %d, %d\n", i,j,k,nx,ny,nz);
			printf( "read_physical_vector(): inconsistent dimensions" );
		}
		
		
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					fscanf( fp, "%lf\n", &w_big[index]);
                    
					index = index + 1;
                }
			}
		}
        
        fclose( fp );
    }
    
    
    
    MPI_Scatterv( u_big,  par->recvcounts_u, par->displs_u, MPI_DOUBLE,
                 u_phys, par->local_size_u,             MPI_DOUBLE,
                 0, MPI_COMM_WORLD );
    MPI_Scatterv( v_big,  par->recvcounts_v, par->displs_v, MPI_DOUBLE,
                 v_phys, par->local_size_v,             MPI_DOUBLE,
                 0, MPI_COMM_WORLD );
    MPI_Scatterv( w_big,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
                 w_phys, par->local_size_wp,             MPI_DOUBLE,
                 0, MPI_COMM_WORLD );
    
    if ( par->my_rank == 0 ) {
        fftw_free( u_big );
        fftw_free( v_big );
        fftw_free( w_big );
    }
    printf("Done reading the velocity IC's\n");
}






extern void bc_dump_physical_vector(
									double *bc_u_phys,
									double *bc_v_phys,
									double *bc_w_phys,
									const char *basename,
									parameters *par )
{
    int nx, ny, nz, total_local_size, total_size;
    double *bc_u_big, *bc_v_big, *bc_w_big;
	
	nx = par->nx;
    ny = par->ny;
    nz = par->nz;
	
    total_local_size = par->local_size_tb;
    total_size       = par->local_size_tb*par->nproc;
	
    if ( par->my_rank == 0 ) {
		
        bc_u_big = ( double * )fftw_malloc( sizeof( double )*total_size );
        bc_v_big = ( double * )fftw_malloc( sizeof( double )*total_size );
        bc_w_big = ( double * )fftw_malloc( sizeof( double )*total_size );
		
        
    }
    else {
        bc_u_big = NULL;
        bc_v_big = NULL;
        bc_w_big = NULL;
		
    }
	
    MPI_Gatherv( bc_u_phys, total_local_size,         MPI_DOUBLE,
				bc_u_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
    MPI_Gatherv( bc_v_phys, total_local_size,         MPI_DOUBLE,
				bc_v_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
    MPI_Gatherv( bc_w_phys, total_local_size,         MPI_DOUBLE,
				bc_w_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
    if ( par->my_rank == 0 ) {
		
        FILE *fp;
        char filename[128];
        int i, k, ijk;
		
        sprintf( filename, "./outputdir/%s_bc_vels_dump_it%d.dat",
				basename,
				par->it );
        fp = fopen( filename, "w" );
        if ( fp == NULL ) {
            
			printf( "bc_dump_physical_vector(): cannot write to file" );
		}
		
        fprintf( fp, "DIMENSIONS %d %d %d\n", nx, 2, nz );
		ijk = 0;
		for (i=0; i<=(par->local_nx+1)*par->nproc-1; ++i){
			for (k=0; k<=nz-1; ++k){
				fprintf( fp, "%+e %+e %+e\n",
                        bc_u_big[ijk], bc_v_big[ijk], bc_w_big[ijk] );
				ijk += 1;
			}
		}
		
        fclose( fp );
    }
	
    if ( par->my_rank == 0 ) {
        fftw_free( bc_u_big );
        fftw_free( bc_v_big );
        fftw_free( bc_w_big );
		
    }
}



extern void bc_u_tau_dump_physical_vector(
                                          double *bc_u_tau,
                                          const char *basename,
                                          parameters *par )
{
    int nx, ny, nz, total_local_size, total_size;
    double *bc_u_tau_big;
    
    nx = par->nx;
    ny = par->ny;
    nz = par->nz;
    
    
    
    total_local_size = par->local_size_tb;
    total_size       = total_local_size*par->nproc;
    if ( par->my_rank == 0 ) {
        bc_u_tau_big = ( double * )fftw_malloc( sizeof( double )*total_size );
    }
    else {
        bc_u_tau_big = NULL;
    }
    
    MPI_Gatherv( bc_u_tau, total_local_size,         MPI_DOUBLE,
                bc_u_tau_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
	printf("Done data exchange %d %d %d %d\n", total_local_size-nz, par->local_nx*nz,par->local_nx, nz);
    if ( par->my_rank == 0 ) {
        FILE *fp;
        char filename[128];
        int i, k, ijk,i_n;
        int i_counter;
        i_counter =0;
        sprintf( filename, "./outputdir/%s_bc_u_tau_dump_it%d.dat",
                     basename,
                par->it );
        fp = fopen( filename, "w" );
        if ( fp == NULL ) {
            
            printf( "bc_u_tau_physical_vector(): cannot write to file\n" );
        }
        
        fprintf( fp, "x z u_tau\n");
        ijk = 0;
        
        for (i_n=0; i_n < par->nproc; i_n++){
            
            for (i=0; i<= par->local_nx; ++i){
				
                for (k=0; k<=nz-1; ++k){
					if (i < par->local_nx){
						printf("%d %d %+e %+e %+e\n",i,k, (double) i_counter*par->dx, ((double) k+0.5)*par->dz, bc_u_tau_big[ijk]);
                        fprintf( fp, "%+e %+e %+e\n", (double) i_counter*par->dx, ((double) k+0.5)*par->dz, bc_u_tau_big[ijk]);
                        
					}
					ijk ++;
                    
				}
				if (i< par->local_nx) {
					i_counter ++;
				}
                
            }
        }
        
        fclose( fp );
    }
    
    if ( par->my_rank == 0 ) {
        fftw_free( bc_u_tau_big );
    }
}


extern void bc_read_physical_vector(
									double *bc_u_phys,
									double *bc_v_phys,
									double *bc_w_phys,
									const char *basename,
									parameters *par )
{
    int nx, ny, nz, total_local_size, total_size;
    double *bc_u_big, *bc_v_big, *bc_w_big;
	
	nx = par->nx;
    ny = par->ny;
    nz = par->nz;
	
    total_local_size = par->local_size_tb;
    total_size       = par->local_size_tb*par->nproc;
	
    if ( par->my_rank == 0 ) {
		
        bc_u_big = ( double * )fftw_malloc( sizeof( double )*total_size );
        bc_v_big = ( double * )fftw_malloc( sizeof( double )*total_size );
        bc_w_big = ( double * )fftw_malloc( sizeof( double )*total_size );
        
    }
    else {
        bc_u_big = NULL;
        bc_v_big = NULL;
        bc_w_big = NULL;
		
    }
	
    if ( par->my_rank == 0 ) {
		
        FILE *fp;
        char filename[128];
        int i, j, k, ijk;
		printf("./input/re%d/%s_bc_vels_dump_%d%d%d.dat\n", (int)par->re_start, basename, nx, ny, nz);
		sprintf( filename, "./input/re%d/%s_bc_vels_dump_%d%d%d.dat",
				(int)par->re_start, basename, nx, ny, nz);
        fp = fopen( filename, "r" );
        if ( fp == NULL ){
			printf( "bc_read_physical_vector(): cannot read from file" );
		}
		
        fscanf( fp, "DIMENSIONS %d %d %d\n", &i, &j, &k );
        if ( !( i == par->nx && j == 2 && k == nz ) ){
            printf("i=%d, j=%d, k=%d, but expecting %d, %d, %d\n", i,j,k,par->nx, 2, par->nz);
            printf( "bc_read_physical_vector(): inconsistent dimensions\n" );
		}
		
		ijk = 0;
		for (i=0; i<=(par->local_nx+1)*par->nproc-1; ++i){
			for (k=0; k<=nz-1; ++k){
				fscanf( fp, "%lf %lf %lf\n",
					   &bc_u_big[ijk], &bc_v_big[ijk], &bc_w_big[ijk] );
				ijk += 1;
            }
		}
		
        fclose( fp );
    }
	
	
    MPI_Scatterv( bc_u_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
				 bc_u_phys, total_local_size,         MPI_DOUBLE,
				 0, MPI_COMM_WORLD );
    MPI_Scatterv( bc_v_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
				 bc_v_phys, total_local_size,         MPI_DOUBLE,
				 0, MPI_COMM_WORLD );
    MPI_Scatterv( bc_w_big,  par->recvcounts_tb, par->displs_tb, MPI_DOUBLE,
				 bc_w_phys, total_local_size,         MPI_DOUBLE,
				 0, MPI_COMM_WORLD );
	
    if ( par->my_rank == 0 ) {
        fftw_free( bc_u_big );
        fftw_free( bc_v_big );
        fftw_free( bc_w_big );
		
    }
}





void check_divergence(
                      variables *var,
                      parameters *par,
                      fftwplans *planptr)
{
	
	
	neighbors sh, *shared;
	shared = &sh;
	interpolated interp, *inter;
	inter = &interp;
	
	int i;
	int local_size_ued = par->local_size_ued;
	int local_size_ved = par->local_size_ved;
	
	
	get_neignboring_data(var, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
    
	get_periodic_data(var, shared, par);
	MPI_Barrier(MPI_COMM_WORLD);
	
	inter->u_ued = FFTW_MALLOC(local_size_ued);
	inter->v_ved = FFTW_MALLOC(local_size_ved);
	
	
	if (par->fd_order == 4){
		bigger_array_u_ued__4(1.0, var, shared, inter->u_ued, par);
		bigger_array_v_ved__4(1.0, var, shared, inter->v_ved, par);
		MPI_Barrier(MPI_COMM_WORLD);
		get_divergence__4(inter->u_ued, inter->v_ved, var->w, var->Du, par);
	}
	if (par->fd_order == 2){
		bigger_array_u_ued__2(1.0, var, shared, inter->u_ued, par);
		bigger_array_v_ved__2(1.0, var, shared, inter->v_ved, par);
		MPI_Barrier(MPI_COMM_WORLD);
		get_divergence__2(inter->u_ued, inter->v_ved, var->w, var->Du, par);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	finalize_neignboring(shared, par);
	fftw_free(inter->u_ued);
	fftw_free(inter->v_ved);
	
	
	fftw_execute_r2r(planptr->p1d_invz_wp, var->Du, var->Du);
    
	MPI_Barrier(MPI_COMM_WORLD);
	double Du_local = 0.0;
	double Du_max = 0.0;
	for (i=0; i<=par->local_size_wp-1; ++i){
		if (Du_local <= pow(var->Du[i], 2.0)){
			Du_local = pow(var->Du[i], 2.0);
		}
        
	}
    
	MPI_Reduce( &Du_local, &Du_max, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	par->Du_max = Du_max/(double)par->nproc;
    
	
	MPI_Barrier(MPI_COMM_WORLD);
    
	return;
}





void initialize_flow_field(
                           variables *var0,
                           variables *var1,
                           parameters *par,
                           fftwplans *planptr)
{
	
	neighbors sh, *shared;
	shared = &sh;
	interpolated interp, *inter;
	inter = &interp;
	int i, k, index;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx_start = par->local_nx_start;
	const int local_size_cn = par->local_size_cn;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	int icont=par->icont;
    
	//*********************** reading data and fftw *******************
    if (icont == 0 ) {
	initialize_velocity_field(var0->u_phys, var0->v_phys, var0->w_phys, par);			
	}
	if (icont == 1)  {
	read_physical_vector(var0->u_phys, var0->v_phys, var0->w_phys, par);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// for checking purpose
	par->it = 17;
	dump_physical_vector(var0->u_phys, var0->v_phys, var0->w_phys, par);
if(par->my_rank == 0){
// printf("dump_physical_vector  succeed %d\n",par->it);
}  
	fftw_execute_r2r(planptr->p1d_z_u, var0->u_phys, var0->u);
	for (i=0; i<=par->local_size_u-1; ++i){
		var0->u[i] = var0->u[i]/(double)nz;
	}
    
	fftw_execute_r2r(planptr->p1d_z_v, var0->v_phys, var0->v);
	for (i=0; i<=par->local_size_v-1; ++i){
		var0->v[i] = var0->v[i]/(double)nz;
	}
    
	fftw_execute_r2r(planptr->p1d_z_wp, var0->w_phys, var0->w);
	for (i=0; i<=par->local_size_wp-1; ++i){
		var0->w[i] = var0->w[i]/(double)nz;
	}
if(par->my_rank == 0){
//printf("fft x y kz  succeed %d\n",par->it);
}  	
	//************************ reading data and fftw *******************
	
	
	if (par->les == 0){
		for (i=0; i<=local_size_cn-1; ++i){
			var0->Txx[i] = 0.0;
			var0->Tyy[i] = 0.0;
			var0->Tzz[i] = 0.0;
            
			var0->Txy[i] = 0.0;
			var0->Tyz[i] = 0.0;
			var0->Tzx[i] = 0.0;
		}
	}
	
	
	if (par->bc == 1){
		
		
		//************************ reading boundary data and fftw *******************
		if (par->bc_data == 1){
			bc_read_physical_vector(var0->bc_u_t, var0->bc_v_t, var0->bc_w_t, "top", par);
			bc_read_physical_vector(var0->bc_u_b, var0->bc_v_b, var0->bc_w_b, "bot", par);
			
			// for checking purpose
			par->it = 17;
			bc_dump_physical_vector(var0->bc_u_t, var0->bc_v_t, var0->bc_w_t, "top", par);
			bc_dump_physical_vector(var0->bc_u_b, var0->bc_v_b, var0->bc_w_b, "bot", par);
			
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_u_t, var0->bc_u_t);
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_u_b, var0->bc_u_b);
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_v_t, var0->bc_v_t);
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_v_b, var0->bc_v_b);
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_w_t, var0->bc_w_t);
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_w_b, var0->bc_w_b);
			for (i=0; i<=par->local_size_tb-1; ++i){
				var0->bc_u_t[i] = var0->bc_u_t[i]/(double)nz;
				var0->bc_u_b[i] = var0->bc_u_b[i]/(double)nz;
				var0->bc_v_t[i] = var0->bc_v_t[i]/(double)nz;
				var0->bc_v_b[i] = var0->bc_v_b[i]/(double)nz;
				var0->bc_w_t[i] = var0->bc_w_t[i]/(double)nz;
				var0->bc_w_b[i] = var0->bc_w_b[i]/(double)nz;
			}
		}
		//************************ reading boundary data and fftw *******************
if(par->my_rank == 0){
//printf("get_neignboring_data begin  succeed %d\n",par->it);
}          
		
		get_neignboring_data(var0, shared, par);
		MPI_Barrier(MPI_COMM_WORLD);
if(par->my_rank == 0){
//printf("get_periodic_data begin  succeed %d\n",par->it);
}		
		get_periodic_data(var0, shared, par);
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		init_interpolation(inter, par);
if(par->my_rank == 0){
//printf("init_interpolation  succeed %d\n",par->it);
}	
		if (par->fd_order == 2){
			interpolation__2(1.0, var0, shared, inter, par, planptr);
		}
		if (par->fd_order == 4){
			interpolation__4(1.0, var0, shared, inter, par, planptr);
		}
if(par->my_rank == 0){
//printf("interpolation__4  succeed %d\n",par->it);
}		
		
		if (par->fd_order == 4){
			get_velgrad_tensor_cn(var0, inter, par, planptr);
		}
        
        	if (par->fd_order == 2){
            		get_velgrad_tensor_cn__2(var0, inter, par, planptr);
       	}
/*if(par->my_rank == 0){
printf("get_velgrad_tensor_cn  succeed %d\n",par->it);
}	 */       
		les_channel(var0, inter, par, planptr);
/*if(par->my_rank == 0){
printf("les_channel  succeed %d\n",par->it);
}	*/		
		finalize_neignboring(shared, par);
		finalize_interpolation(inter, par);
		
		
		double *bc_u_t;
		double *bc_u_b;
		bc_u_t = dvector(0, (local_nx+1)*nz-1);
		bc_u_b = dvector(0, (local_nx+1)*nz-1);
		
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
				bc_u_t[index] = var0->u_phys[par->index_u[i][ny-1][k]];
				bc_u_b[index] = var0->u_phys[par->index_u[i][0][k]];
			}
		}
		
		bc__get_bc_dudy(bc_u_t, var0->bc_kappa_t, var0->bc_dudy_t, "top", par);
		bc__get_bc_dudy(bc_u_b, var0->bc_kappa_b, var0->bc_dudy_b, "bot", par);
	
		
		free_dvector(bc_u_b, 0, (local_nx+1)*nz-1);
		free_dvector(bc_u_t, 0, (local_nx+1)*nz-1);
        
		if (par->bc_data == 0){
			// if no boundary data is available, calculate from interior velocity data
            
			bc__get_bc_vel(var0->bc_dudy_t, var0->bc_kappa_t, var0->bc_kappa_t,
						   var0->bc_u_t, var0->bc_v_t, var0->bc_w_t, "top", par);
			bc__get_bc_vel(var0->bc_dudy_b, var0->bc_kappa_b, var0->bc_kappa_b,
						   var0->bc_u_b, var0->bc_v_b, var0->bc_w_b, "bot", par);
            
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_u_t, var0->bc_u_t);
			fftw_execute_r2r(planptr->p1d_z_bc_tb, var0->bc_u_b, var0->bc_u_b);
			for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
				var0->bc_u_t[i] = var0->bc_u_t[i]/(double)(nz);
				var0->bc_u_b[i] = var0->bc_u_b[i]/(double)(nz);
			}
			

		}
        
		
		
	}
}




void initialize_bc_flow_field(
                              variables *var,
                              parameters *par,
                              fftwplans *planptr)
{
	
	int i, k;
	const int local_nx = par->local_nx;
	const int nz = par->nz;
	const int local_nx_start = par->local_nx_start;
	const int start = local_nx_start;
	const int end = local_nx_start + local_nx;
	
	const double dx = par->dx;
	const double dz = par->dz;
	
	double x, z;
    
	for (i=start; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			var->bc_v_t[par->index_tb[i][k]] = 0.0;
			var->bc_w_t[par->index_tb[i][k]] = 0.0;
		}
	}
	
	for (i=start; i<=end; ++i){
		for (k=0; k<=nz-1; ++k){
			x = dx*i;
			z = dz*k;
			var->bc_u_t[par->index_tb[i][k]] = 100.0*cos(PI*x)*sin(2*PI*z);
			
			var->bc_u_b[par->index_tb[i][k]] = - 100.0*cos(PI*x)*sin(2*PI*z);
		}
	}
	fftw_execute_r2r(planptr->p1d_z_bc_tb, var->bc_u_t, var->bc_u_t);
	fftw_execute_r2r(planptr->p1d_z_bc_tb, var->bc_v_t, var->bc_v_t);
	fftw_execute_r2r(planptr->p1d_z_bc_tb, var->bc_w_t, var->bc_w_t);
	fftw_execute_r2r(planptr->p1d_z_bc_tb, var->bc_u_b, var->bc_u_b);
	fftw_execute_r2r(planptr->p1d_z_bc_tb, var->bc_v_b, var->bc_v_b);
	fftw_execute_r2r(planptr->p1d_z_bc_tb, var->bc_w_b, var->bc_w_b);
	for (i=0; i<=(nz)*(local_nx+1)-1; ++i){
		var->bc_u_t[i] = var->bc_u_t[i]/(double)(nz);
		var->bc_v_t[i] = var->bc_v_t[i]/(double)(nz);
		var->bc_w_t[i] = var->bc_w_t[i]/(double)(nz);
		var->bc_u_b[i] = var->bc_u_b[i]/(double)(nz);
		var->bc_v_b[i] = var->bc_v_b[i]/(double)(nz);
		var->bc_w_b[i] = var->bc_w_b[i]/(double)(nz);
	}
}




extern void get_mass_flux(
                          double *u,
                          parameters *par )
{
	
	int i, j;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	
	double local_u_y_mean;
	double u_y_mean;
    
	
	local_u_y_mean = 0.0;
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			local_u_y_mean += u[par->index_u[i][j][0]];
		}
	}
	local_u_y_mean = local_u_y_mean/(double)(local_nx*ny);
	
	MPI_Reduce(&local_u_y_mean, &u_y_mean,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	
	if (my_rank == 0){
		par->mass_flux_start = u_y_mean/(double)par->nproc;
		printf("Um = %e (%e)\n", par->mass_flux_start, par->mass_flux_aim);
	}
	
	MPI_Bcast( &par->mass_flux_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
}

extern void dump_plane_data_plot(
                                 double *plane,
                                 const char *basename,
                                 parameters *par)
{
	
	int total_local_size, total_size;
    double *big;
    const int ny = par->ny;
	
    total_local_size = par->local_nx*ny;
    total_size       = par->local_nx*ny*par->nproc;
	
    if ( par->my_rank == 0 ) {
        big = ( double * )fftw_malloc( sizeof( double )*total_size );
    }
    else {
        big = NULL;
    }
    
    
    MPI_Gatherv(plane, total_local_size,         MPI_DOUBLE,
				big,  par->recvcounts_mean_2d, par->displs_mean_2d, MPI_DOUBLE,
				0, MPI_COMM_WORLD);
	
    if ( par->my_rank == 0 ) {
		
        FILE *fp;
        char filename[128];
        int i, j, index;
		double x, y;
		
        sprintf( filename, "./outputdir/%s_it%d.txt",
				basename,
				par->it );
        fp = fopen( filename, "w" );
        if ( fp == NULL ) {
			printf( "dump_mean_vector(): cannot write to file" );
		}
		
        fprintf( fp, "#x	y\n");
		
		index = 0;
		for (i=0; i<=par->local_nx*par->nproc-1; ++i){
			for (j=0; j<=ny-1; ++j){
				x = (i+0.5)*par->dx;
				y = (j+0.5)*par->dy;
				fprintf( fp, "%+e	%+e	%+e\n",
                        x, y, big[index]);
				index += 1;
			}
			fprintf( fp, "\n");
		}
        fclose( fp );
    }
	
    if ( par->my_rank == 0 ) {
        fftw_free( big );
    }
	
	return;
}


extern void dump_physical_vector_plot(
									  double *u_phys,
									  double *v_phys,
									  double *w_phys,
									  parameters *par )
{
	int i, j, index;
	const int ny = par->ny;
	const int local_nx = par->local_nx;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *u_k0, *v_k0, *w_k0;
	
	u_k0 = FFTW_MALLOC(ny*local_nx);
	v_k0 = FFTW_MALLOC(ny*local_nx);
	w_k0 = FFTW_MALLOC(ny*local_nx);
	
	index = 0;
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			u_k0[index] = u_phys[par->index_u[i][j][0]];
			
			if (j==0){
				v_k0[index] = 0.0;
			}else{
				v_k0[index] = v_phys[par->index_v[i][j][0]];
			}
			w_k0[index] = w_phys[par->index_wp[i][j][0]];
			index += 1;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	dump_plane_data_plot(u_k0, "u_field", par);
	dump_plane_data_plot(v_k0, "v_field", par);
	dump_plane_data_plot(w_k0, "w_field", par);
	
	fftw_free(u_k0);
	fftw_free(v_k0);
	fftw_free(w_k0);
	
	return;
}

// generate initial random velocity field by code2016.11.10
extern void initialize_velocity_field(  
                                 double *u_phys,
                                 double *v_phys,
                                 double *w_phys,
                                 parameters *par )
{
	int nx,ny,nz;
	int total_size_u;
	int total_size_v;
	int total_size_wp;
	double *u_big,*v_big,*w_big;

	nx=par->nx;
	ny=par->ny;
	nz=par->nz;

	total_size_u=nx*ny*nz;
	total_size_v=nx*(ny-1)*nz;
	total_size_wp=nx*ny*nz;
	printf("initial random field with laminar flow");
	    if ( par->my_rank == 0 ) {
        u_big = (double *)fftw_malloc(sizeof(double)*total_size_u);
        v_big = (double *)fftw_malloc(sizeof(double)*total_size_v);
        w_big = (double *)fftw_malloc(sizeof(double)*total_size_wp);
    }
    else {
        u_big = NULL;
        v_big = NULL;
        w_big = NULL;
    }

	if( par->my_rank == 0 ){
	    int i,j,k,index,number;
		double y;
// u 11.10 wan shang ji xu gao ding
		//random seed.
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
					y=(j+0.5-ny/2.0)*par->dy;
				for (k=0; k<=nz-1; ++k){
					number=rand()%101;
					if(i==0 && k==0){
						printf("number=%f\n",(number/100.0-0.5)/10.0);
					}
					u_big[index]=(1.0-y*y)*3.0*par->mass_flux_aim*\
								(1+(number/100.0-0.5)/10.0);                
					index = index + 1;
                }
			}
		}
// v
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=1; j<=ny-1; ++j){
				y=(j+0.5-ny/2.0)*par->dy;
				for (k=0; k<=nz-1; ++k){
					v_big[index]=0.0;
					number=rand()%101;
					v_big[index]=(1.0-y*y)*3.0*par->mass_flux_aim*\
								(0.0+(number/100.0-0.5)/10.0); 
            		index = index + 1;
                }
			}
		}
        
// w
		index = 0;
		for (i=0; i<=nx-1; ++i){
			for (j=0; j<=ny-1; ++j){
				for (k=0; k<=nz-1; ++k){
					w_big[index]=0.0; 
					number=rand()%101;
					w_big[index]=(1.0-y*y)*3.0*par->mass_flux_aim*\
								(0.0+(number/100.0-0.5)/10.0); 
					index = index + 1;
                }
			}
		}
	}

	MPI_Scatterv( u_big,  par->recvcounts_u, par->displs_u, MPI_DOUBLE,
                 u_phys, par->local_size_u,             MPI_DOUBLE,
                 0, MPI_COMM_WORLD );
    MPI_Scatterv( v_big,  par->recvcounts_v, par->displs_v, MPI_DOUBLE,
                 v_phys, par->local_size_v,             MPI_DOUBLE,
                 0, MPI_COMM_WORLD );
    MPI_Scatterv( w_big,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
                 w_phys, par->local_size_wp,             MPI_DOUBLE,
                 0, MPI_COMM_WORLD );
    
    if ( par->my_rank == 0 ) {
        fftw_free( u_big );
        fftw_free( v_big );
        fftw_free( w_big );
    }
    printf("Done generating initial velocity field\n");
}
