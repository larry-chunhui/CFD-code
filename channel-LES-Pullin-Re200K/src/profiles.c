/* profiles.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 */

#include <math.h>
#include <mpi.h>
#include "definitions.h"
#include "nrutil.h"
#include "profiles.h"


extern void skin_friction(
                          double *u,
                          parameters *par )
{
    
	int i;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double nu = 1/par->re;
	
	double local_mean_dudy;
	double mean_dudy;
	double dudy_t;
	double dudy_b;
	
	dudy_t = 0.0;
	dudy_b = 0.0;
	for (i=start; i<=end-1; ++i){
        
		if (par->fd_order == 2){
			dudy_t += (0 - 9*u[par->index_u[i][ny-1][0]] + u[par->index_u[i][ny-2][0]])/(3.0*par->dy);
			dudy_b += (0 + 9*u[par->index_u[i][0][0]] - u[par->index_u[i][1][0]])/(3.0*par->dy);
		}
		if (par->fd_order == 4){
			dudy_t += (0 - 225*u[par->index_u[i][ny-1][0]] + 50*u[par->index_u[i][ny-2][0]] - 9*u[par->index_u[i][ny-3][0]] )/(60.0*par->dy);
			dudy_b += (0 + 225*u[par->index_u[i][0][0]] - 50*u[par->index_u[i][1][0]] + 9*u[par->index_u[i][2][0]])/(60.0*par->dy);
		}
		
	}
	
	// take a local mean of dudy
	local_mean_dudy = (dudy_b - dudy_t)/(double)(2*local_nx);
	
	MPI_Reduce(&local_mean_dudy, &mean_dudy,
               1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	if (my_rank == 0){
		mean_dudy = mean_dudy/(double)par->nproc;
		// friction velocity
		par->u_tau = sqrt(fabs(nu*mean_dudy));
	}
	
	MPI_Bcast( &par->u_tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
}


extern void dump_bc_statistics(
                               variables *var,
                               parameters *par)
{
	
	FILE *fp;
	char filename[100];
	
	int i, k, index;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	
	if (my_rank == 0 && par->it%par->stat == 0){
		

		sprintf( filename, "./%s/bc_profile_%d.dat", par->outputdir, par->it);
		fp = fopen( filename, "a" );
		
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				index = par->index_tb[i][k];
				fprintf( fp, "%d	%d	%+e	%+e	%+e	%+e	%+e	%+e\n",
						i,
						k,
						var->bc_kappa_t[index],
						var->bc_dudy_t[index],
						var->bc_txy_t[index],
						var->bc_kappa_b[index],
						var->bc_dudy_b[index],
						var->bc_txy_b[index]);
			}
			fprintf( fp, "\n");
		}
		
		fclose( fp );
	}
	
}




extern void skin_friction__bc(
                              variables *var,
                              parameters *par )
{
    
	int i, k, index;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double nu = 1/par->re;
	
	double local_mean_dudy;
	double local_mean_dudy_vw_t;
	double local_mean_dudy_vw_b;
	double local_mean_kappa;
	double local_mean_bcK_t;
	double local_mean_bctxy_t;
	double local_mean_bcK_b;
	double local_mean_bctxy_b;
	double mean_dudy;
	double mean_dudy_vw_t;
	double mean_dudy_vw_b;
	double mean_kappa;
	double mean_bcK_t;
	double mean_bctxy_t;
	double mean_bcK_b;
	double mean_bctxy_b;
	double dudy_t;
	double dudy_b;
	
	dudy_t = 0.0;
	dudy_b = 0.0;
	local_mean_dudy_vw_t = 0.0;
	local_mean_dudy_vw_b = 0.0;
	local_mean_kappa = 0.0;
	local_mean_bcK_t = 0.0;
	local_mean_bcK_b = 0.0;
	local_mean_bctxy_t = 0.0;
	local_mean_bctxy_b = 0.0;
    
	for (i=start; i<=end-1; ++i){
		for (k=0; k<=nz-1; ++k){
			index = par->index_tb[i][k];
			dudy_b += var->bc_dudy_b[index];
			dudy_t += var->bc_dudy_t[index];
			local_mean_kappa += var->bc_kappa_t[index];
			local_mean_kappa += var->bc_kappa_b[index];
			local_mean_bcK_t += var->bc_K_t[index];
			local_mean_bcK_b += var->bc_K_b[index];
			local_mean_bctxy_t += var->bc_txy_t[index];
			local_mean_bctxy_b += var->bc_txy_b[index];
			
			// at virtual wall
			index = par->index_cn[i][0][k];
			local_mean_dudy_vw_b += var->dudy[index];
			index = par->index_cn[i][ny][k];
			local_mean_dudy_vw_t += var->dudy[index];
		}
	}
	
	// take a local mean of dudy
	local_mean_dudy = (dudy_b - dudy_t)/(double)(2*local_nx*nz);
	local_mean_dudy_vw_t = local_mean_dudy_vw_t/(double)(local_nx*nz);
	local_mean_dudy_vw_b = local_mean_dudy_vw_b/(double)(local_nx*nz);
	local_mean_kappa = local_mean_kappa/(double)(2*local_nx*nz);
	local_mean_bcK_t = local_mean_bcK_t/(double)(local_nx*nz);
	local_mean_bctxy_t = local_mean_bctxy_t/(double)(local_nx*nz);
	local_mean_bcK_b = local_mean_bcK_b/(double)(local_nx*nz);
	local_mean_bctxy_b = local_mean_bctxy_b/(double)(local_nx*nz);
	
	MPI_Reduce(&local_mean_dudy, &mean_dudy,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_dudy_vw_t, &mean_dudy_vw_t,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_dudy_vw_b, &mean_dudy_vw_b,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_kappa, &mean_kappa,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_bcK_t, &mean_bcK_t,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_bctxy_t, &mean_bctxy_t,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_bcK_b, &mean_bcK_b,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	MPI_Reduce(&local_mean_bctxy_b, &mean_bctxy_b,
			   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	if (my_rank == 0){
		mean_dudy = mean_dudy/(double)par->nproc;
		
		mean_dudy_vw_t = mean_dudy_vw_t/(double)par->nproc;
		mean_dudy_vw_b = mean_dudy_vw_b/(double)par->nproc;
		
		par->mean_kappa = mean_kappa/(double)par->nproc;
		
		par->mean_bcK_t = mean_bcK_t/(double)par->nproc;
		par->mean_bctxy_t = mean_bctxy_t/(double)par->nproc;
		
		par->mean_bcK_b = mean_bcK_b/(double)par->nproc;
		par->mean_bctxy_b = mean_bctxy_b/(double)par->nproc;
		
		// friction velocity
		par->u_tau = sqrt(fabs(nu*mean_dudy));
		
		par->mean_dudy = mean_dudy;
		par->expected_dpdx_t = nu*mean_dudy_vw_t - par->mean_bctxy_t;
		par->expected_dpdx_b = nu*mean_dudy_vw_b - par->mean_bctxy_b;
	}
	
	MPI_Bcast( &par->u_tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->mean_kappa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->mean_bcK_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->mean_bctxy_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->mean_bcK_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->mean_bctxy_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->expected_dpdx_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &par->expected_dpdx_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	
}


//10.2添加w_mean, u_rms,v_rms,w_rms,<uv>
extern void statistics_u(
                         double *u,
                         double *v,
			    double *w,
                         parameters *par )
{
	
	int i, j;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double delta = par->Ly/2.0;
	const double nu = 1/par->re, u_tau = par->u_tau;
	
	FILE *fp;
	char filename[100];
	
	double local_u_y_mean, local_u_y_centerline;
	double u_y_mean, u_y_centerline;

//-------------------------------	
// add on 2016.10.2
	double local_u_mean, scratch;
	double local_v_mean;
	double local_w_mean;

	double *u_mean = NULL;
	double *v_mean = NULL;
	double *w_mean = NULL;

	double local_u_rms,local_v_rms,local_w_rms,local_uv_mean;
	double *u_rms = NULL;
	double *v_rms = NULL;
	double *w_rms = NULL;
	double *uv_mean = NULL;
	double ek_t;

		u_mean = dvector(0, ny-1);
		v_mean = dvector(0, ny-1);
		w_mean = dvector(0, ny-1);

//-------------------------------	
// add on 2016.10.2
//-------------------------------	
	if (my_rank == 0){
		u_rms = dvector(0,ny-1);
		v_rms = dvector(0,ny-1);
		w_rms = dvector(0,ny-1);
		uv_mean=dvector(1,ny-1);
//printf("allocate uvw mean and rms succeed\n");


	}
//u（y）i have not change this one.	
	for (j=0; j<=ny-1; ++j){
		local_u_mean = 0.0;
		for (i=start; i<=end-1; ++i){
			local_u_mean += u[par->index_u[i][j][0]];
		}
		local_u_mean = local_u_mean/(double)(local_nx);
		MPI_Reduce(&local_u_mean, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			u_mean[j] = scratch/(double)par->nproc;
		}
	}

	for (j=1; j<=ny-1; ++j){
		local_v_mean = 0.0;
		for (i=start; i<=end-1; ++i){
			local_v_mean += v[par->index_v[i][j][0]];
		}
		local_v_mean = local_v_mean/(double)(local_nx);
		MPI_Reduce(&local_v_mean, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			v_mean[j] = scratch/(double)par->nproc;
		}
		
	}
		if (my_rank == 0){
			v_mean[0] = 0.00;
		}

	for (j=0; j<=ny-1; ++j){
		local_w_mean = 0.0;
		for (i=start; i<=end-1; ++i){
			local_w_mean += w[par->index_wp[i][j][0]];
		}
		local_w_mean = local_w_mean/(double)(local_nx);
		MPI_Reduce(&local_w_mean, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			w_mean[j] = scratch/(double)par->nproc;
		}
		
	}

    	MPI_Bcast(u_mean, ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    	MPI_Bcast(v_mean, ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    	MPI_Bcast(w_mean, ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);      //all succeed
    	MPI_Barrier(MPI_COMM_WORLD);

//		printf("my_rank=%d,u_mean[%d]=%e\n",par->my_rank,ny/2,u_mean[ny/2]);
//		printf("my_rank=%d,v_mean[%d]=%e\n",par->my_rank,ny/2,v_mean[ny/2]);
//		printf("my_rank=%d,w_mean[%d]=%e\n",par->my_rank,ny/2,w_mean[ny/2]);
// u_rms
	for (j=0; j<=ny-1; ++j){
		local_u_rms = 0.0;

		for (i=start; i<=end-1; ++i){
			local_u_rms += pow((u[par->index_u[i][j][0]]-u_mean[j]),2.00);  //一会在改正
		}

		local_u_rms = local_u_rms/(double)(local_nx);

		MPI_Reduce(&local_u_rms, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );       //

		if (my_rank == 0 ){
			u_rms[j] =sqrt( scratch/(double)par->nproc);
		}		
	}
// v_rms
	for (j=1; j<=ny-1; ++j){
		local_v_rms = 0.0;
		for (i=start; i<=end-1; ++i){
			local_v_rms += pow((v[par->index_v[i][j][0]]-v_mean[j]),2.00);  //一会在改正
		}
		local_v_rms = local_v_rms/(double)(local_nx);
		MPI_Reduce(&local_v_rms, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			v_rms[j] = sqrt(scratch/(double)par->nproc);
		}
		
	}
// w_rms
	for (j=0; j<=ny-1; ++j){
		local_w_rms = 0.0;
		for (i=start; i<=end-1; ++i){
			local_w_rms += pow((w[par->index_wp[i][j][0]]-w_mean[j]),2.00);  //一会在改正
		}
		local_w_rms = local_w_rms/(double)(local_nx);
		MPI_Reduce(&local_w_rms, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			w_rms[j] = sqrt( scratch/(double)par->nproc);
		}
		
	}

// <uv>10.2 add
	for (j=1; j<=ny-1; ++j){
		local_uv_mean = 0.0;
		for (i=start; i<=end-1; ++i){
			local_uv_mean +=(u[par->index_u[i][j-1][0]]-u_mean[j-1]+   \
							 u[par->index_u[i][j][0]]-u_mean[j])/2.00* \
							(v[par->index_v[i][j][0]]-v_mean[j]);  //
		}
		local_uv_mean = local_uv_mean/(double)(local_nx);
		MPI_Reduce(&local_uv_mean, &scratch,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			uv_mean[j] = scratch/(double)par->nproc;
		}
		
	}
// 输出全场平均湍动能和雷诺应力uv
		if (my_rank == 0){
			ek_t=0.0;
			for(j=1;j<=ny-1;++j){
				ek_t=ek_t+pow(u_rms[j],2.0)+\
						  pow(v_rms[j],2.0)+\
						  pow(w_rms[j],2.0);
			}
			ek_t=ek_t/(ny-1);
		}
//对u沿xyz三个方向做平均	
	local_u_y_mean = 0.0;
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			local_u_y_mean += u[par->index_u[i][j][0]];
		}
	}
	local_u_y_mean = local_u_y_mean/(double)(local_nx*ny);
	
	MPI_Reduce(&local_u_y_mean, &u_y_mean,
               1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    
	local_u_y_centerline = 0.0;
	for (i=start; i<=end-1; ++i){
		if (ny%2 == 0){
			local_u_y_centerline += (u[par->index_u[i][ny/2-1][0]] + u[par->index_u[i][ny/2][0]])/2.0;
		}else{
			local_u_y_centerline += u[par->index_u[i][ny/2][0]];
		}
	}
	local_u_y_centerline = local_u_y_centerline/(double)local_nx;
	
	MPI_Reduce(&local_u_y_centerline, &u_y_centerline,
               1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	
	if (my_rank == 0){
		u_y_mean = u_y_mean/(double)par->nproc;
		u_y_centerline = u_y_centerline/(double)par->nproc;
		par->mass_flux_curr = u_y_mean;
		printf("Um = %e (%e)	Uc = %e\n", u_y_mean, par->mass_flux, u_y_centerline);
	}
	
	MPI_Bcast( &par->mass_flux_curr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	
	
	if (my_rank == 0 && par->it%par->stat == 0){
// 此处计算的仅仅是解析的雷诺应力部分		
//		sprintf( filename, "./outputdir/u_y_%d.dat",par->it);
//		fp = fopen( filename, "w" );
//		fprintf( fp, "variables=y,u_mean,v_mean,w_mean,u_rms,v_rms,w_rms,uv_mean\n");       
//		for (j=1; j<=ny/2; ++j){
//			fprintf( fp, "%e,%e,%e,%e,%e,%e,%e,%e\n",
//                  par->dy*(j),
//                    (u_mean[j-1]+u_mean[j])/2.0/par->u_tau,
//					 v_mean[j]/par->u_tau,
//					(w_mean[j-1]+w_mean[j])/2.0/par->u_tau,
//					(u_rms[j-1]+u_rms[j])/2.0/par->u_tau,
//					 v_rms[j]/par->u_tau,
//					(w_rms[j-1]+w_rms[j])/2.0/par->u_tau,
//					 uv_mean[j])/par->u_tau;
//		}     
//		fclose( fp );
//---------------------------------------
		sprintf( filename, "./outputdir/u_profile_%d.dat",par->it);
		fp = fopen( filename, "w" );
		fprintf( fp, "variables=y,u_mean\n");       
		for (j=0; j<=ny-1; ++j){
			fprintf( fp, "%e,%e\n",
                    par->dy*(j+0.5),
                    u_mean[j]);
		}     
		fclose( fp );
		
	}
	
	if (my_rank == 0 && par->it%par->choice == 0){
		sprintf( filename, "./outputdir/statistics_u.dat");
		fp = fopen( filename, "a" );   
		if (par->it == 0){  //change on 11.14
			fprintf( fp, "variables=time,ekt,uc,um,u_tau,re_tau,uv2,uv4,K,dudy,dPdx,dpdx_t,dpdx_b,Cf\n");
		}
		
		fprintf( fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",
                par->tott,
				ek_t,
                u_y_centerline,
                u_y_mean,
                u_tau,
                u_tau*delta/nu,
				uv_mean[2],
				uv_mean[4],
                par->mean_kappa,
                par->mean_dudy,
                par->dPdx,
                par->expected_dpdx_t,
                par->expected_dpdx_b,
                2.0*u_tau*u_tau/u_y_centerline/u_y_centerline);
        
		fclose( fp );
	}

		free_dvector(u_mean, 0, ny-1);
		free_dvector(v_mean, 0, ny-1);
		free_dvector(w_mean, 0, ny-1);
	
	if (my_rank == 0){
		free_dvector(u_rms, 0, ny-1);
		free_dvector(v_rms, 0, ny-1);
		free_dvector(w_rms, 0, ny-1);
		free_dvector(uv_mean,1,ny-1);
	}
    
}







extern void mean_dissipation_ratio(
                                   double *res_diss,
                                   double *sub_diss,
                                   parameters *par )
{
	
	int i, j, k;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	
	FILE *fp;
	char filename[100];
	
	double local_res, local_sgs;
	double global_res, global_sgs;
	double *ratio_mean = NULL;
	double *nusijsij = NULL;
	double *sijtij = NULL;
	double nusijsij_avg=0.0;
	double sijtij_avg=0.0;
    
	if (my_rank == 0){
		ratio_mean = dvector(0, ny);
 		nusijsij = dvector(0, ny);
		sijtij = dvector(0, ny);
        
	}
	
	for (j=0; j<=ny; ++j){
		local_res = 0.0;
		global_res = 0.0;
		local_sgs = 0.0;
		global_sgs = 0.0;
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				local_res += res_diss[par->index_cn[i][j][k]];
				local_sgs += sub_diss[par->index_cn[i][j][k]];
				
			}
		}
		
		local_res = local_res/(double)(local_nx*nz);
		local_sgs = local_sgs/(double)(local_nx*nz);
		
		MPI_Reduce(&local_res, &global_res,
				   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce(&local_sgs, &global_sgs,
				   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			global_sgs /= (double) (par->nproc);
			global_res /= (double) (par->nproc);
			ratio_mean[j] = global_sgs/(global_res+global_sgs);
			nusijsij[j] = global_res;
			sijtij[j] = global_sgs;
		}
	}
	
	if (my_rank == 0){
		
		sprintf( filename, "./%s/dissipation_ratio_profile_%d.dat", par->outputdir, par->it);
		fp = fopen( filename, "a" );
		
		for (j=0; j<=ny; ++j){
			fprintf( fp, "%+e %+e\n",
					par->dy*j,
					ratio_mean[j]);
		}
		
		fclose( fp );
        
        for (j=0; j<=ny; ++j){
            nusijsij_avg += nusijsij[j];
            sijtij_avg += sijtij[j];
        }
        nusijsij_avg /= (double) (ny+1);
        sijtij_avg /= (double) (ny+1);
        
		sprintf( filename, "./%s/dissipation_hist.dat", par->outputdir);
        fp = fopen( filename, "a" );
        
		fprintf( fp, "%d     %+e     %+e    %+e    \n", par->it, nusijsij_avg, sijtij_avg, nusijsij_avg+sijtij_avg);
        
        fclose( fp );
	}
	
	if (my_rank == 0){
		free_dvector(ratio_mean, 0, ny);
		free_dvector(nusijsij, 0, ny);
		free_dvector(sijtij, 0, ny);
	}
	
}





extern void mean_profile_eta(
                             double *res_diss,
                             parameters *par )
{
	
	int i, j, k;
	const int my_rank = par->my_rank;
	const int local_nx = par->local_nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + local_nx;
	const double nu = 1.0/par->re;
    const double delta = pow( ( par->dx )
                             *( par->dy )
                             *( par->dz ), 1.0/3.0 );
	
	FILE *fp;
	char filename[100];
	
	double local_eta_mean, scratch;
	double *eta_mean = NULL;
	
	if (my_rank == 0){
		eta_mean = dvector(0, ny);
	}
	
	for (j=0; j<=ny; ++j){
		local_eta_mean = 0.0;
		for (i=start; i<=end-1; ++i){
			for (k=0; k<=nz-1; ++k){
				local_eta_mean 
				+= pow(nu, 3.0/4.0)/pow(res_diss[par->index_cn[i][j][k]], 1.0/4.0);
			}
		}
		local_eta_mean = local_eta_mean/(double)(local_nx*nz);
		MPI_Reduce(&local_eta_mean, &scratch,
				   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		if (my_rank == 0){
			eta_mean[j] = scratch/(double)par->nproc;
		}
	}
    
	
	if (my_rank == 0){
		
		sprintf( filename, "./%s/eta_mean_profile_%d.dat", par->outputdir, par->it);
		fp = fopen( filename, "a" );
		
		for (j=0; j<=ny; ++j){
			fprintf( fp, "%+e %+e\n",
					par->dy*j,
					delta/eta_mean[j]);
		}
		
		fclose( fp );
		
	}
	
	
	if (my_rank == 0){
		free_dvector(eta_mean, 0, ny);
	}
	
}

