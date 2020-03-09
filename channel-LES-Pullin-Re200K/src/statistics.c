/* statistics.c
 *
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Calculates the statistical quantities such as turbulent intensities. 
 */


#include <mpi.h>
#include <math.h>
#include "definitions.h"
#include "nrutil.h"
#include "statistics.h"
#include "interpolation__4.h"
#include "interpolation__2.h"


extern void statistics_les(
                           interpolated *inter,
                           fftwplans *planptr,
                           parameters *par )
{
	
	int	i, j, k;
	
	const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
	
	const int local_size_wp = par->local_size_wp;
	const int	start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	double *scratch;
	double *u_xz;
	double *v_xz;
	double *w_xz;
	double *u_phys;
	double *v_phys;
	double *w_phys;
	
	double *Txx_xz;
	double *Tyy_xz;
	double *Tzz_xz;
	double *Txy_xz;
	
	double *Txx_phys;
	double *Tyy_phys;
	double *Tzz_phys;
	double *Txy_phys;
	
	
	scratch = FFTW_MALLOC(par->local_size_wp);
	
	if (par->my_rank == 0){
		u_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		v_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		w_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		u_phys = FFTW_MALLOC(nx*ny*nz);
		v_phys = FFTW_MALLOC(nx*ny*nz);
		w_phys = FFTW_MALLOC(nx*ny*nz);
		
		Txx_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		Tyy_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		Tzz_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		Txy_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		Txx_phys = FFTW_MALLOC(nx*ny*nz);
		Tyy_phys = FFTW_MALLOC(nx*ny*nz);
		Tzz_phys = FFTW_MALLOC(nx*ny*nz);
		Txy_phys = FFTW_MALLOC(nx*ny*nz);
		
		for (i=0; i<=nx*ny*2*(nz/2+1)-1; ++i){
			u_xz[i] = 0.0;
			v_xz[i] = 0.0;
			w_xz[i] = 0.0;
			Txx_xz[i] = 0.0;
			Tyy_xz[i] = 0.0;
			Tzz_xz[i] = 0.0;
			Txy_xz[i] = 0.0;
		}
		
	}else{
		u_xz = NULL;
		v_xz = NULL;
		w_xz = NULL;
		u_phys =NULL;
		v_phys = NULL;
		w_phys = NULL;
		
		Txx_xz = NULL;
		Tyy_xz = NULL;
		Tzz_xz = NULL;
		Txy_xz = NULL;
		Txx_phys =NULL;
		Tyy_phys = NULL;
		Tzz_phys = NULL;
		Txy_phys = NULL;
	}
	
	// ***************** resolved stress ***************************
	// get physical data defined at center
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (par->fd_order == 4){
					scratch[par->index_wp[i][j][k]] \
					= 9.0/8.0*inter->u_phys_ct1[par->index_ct[i][j][k]] \
					- 1.0/8.0*inter->u_phys_ct3[par->index_ct[i][j][k]];
				}
				if (par->fd_order == 2){
					scratch[par->index_wp[i][j][k]] \
					= inter->u_phys_ct1[par->index_ct[i][j][k]];
				}
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				u_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (par->fd_order == 4){
					scratch[par->index_wp[i][j][k]] \
					= 9.0/8.0*inter->v_phys_ct1[par->index_ct[i][j][k]] \
					- 1.0/8.0*inter->v_phys_ct3[par->index_ct[i][j][k]];
					
				}
				if (par->fd_order == 2){
					scratch[par->index_wp[i][j][k]] \
					= inter->v_phys_ct1[par->index_ct[i][j][k]];
				}
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				v_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				scratch[par->index_wp[i][j][k]] \
				= inter->w_phys_ct[par->index_ct[i][j][k]];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				w_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	// ***************** resolved stress ***************************
	
	
	// ***************** SGS stress ********************************
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				scratch[par->index_wp[i][j][k]] \
				= inter->Txx_ct[par->index_ct[i][j][k]];
			}
		}
	}
	fftw_execute_r2r(planptr->p1d_invz_wp, scratch, scratch);
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				Txx_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				scratch[par->index_wp[i][j][k]] \
				= inter->Tyy_ct[par->index_ct[i][j][k]];
			}
		}
	}
	fftw_execute_r2r(planptr->p1d_invz_wp, scratch, scratch);
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				Tyy_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				scratch[par->index_wp[i][j][k]] \
				= inter->Tzz_ct[par->index_ct[i][j][k]];
			}
		}
	}
	fftw_execute_r2r(planptr->p1d_invz_wp, scratch, scratch);
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				Tzz_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	
	for (i=0; i<=local_size_wp-1; ++i){
		scratch[i] = 0.0;
	}
	
    
    if (par->fd_order == 4){
        interpolate_cn2wp__4(inter->Txy_cn, scratch, par);
    }
    if (par->fd_order == 2){
        interpolate_cn2wp__2(inter->Txy_cn, scratch, par);
    }
    
	fftw_execute_r2r(planptr->p1d_invz_wp, scratch, scratch);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
				Txy_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
				0, MPI_COMM_WORLD );
	
	
	fftw_free(scratch);
	// ***************** SGS stress ********************************
	
	
	if (par->my_rank == 0){
		
		forward_xz(u_phys, u_xz, planptr, par);
		forward_xz(v_phys, v_xz, planptr, par);
		forward_xz(w_phys, w_xz, planptr, par);
		forward_xz(Txx_phys, Txx_xz, planptr, par);
		forward_xz(Tyy_phys, Tyy_xz, planptr, par);
		forward_xz(Tzz_phys, Tzz_xz, planptr, par);
		forward_xz(Txy_phys, Txy_xz, planptr, par);
		
		// get reystress
		dump_reystress_les(u_xz, v_xz, w_xz,
						   Txx_xz, Tyy_xz, Tzz_xz, Txy_xz, par);
		
		
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (par->my_rank == 0){
		fftw_free(u_phys);
		fftw_free(v_phys);
		fftw_free(w_phys);
		fftw_free(u_xz);
		fftw_free(v_xz);
		fftw_free(w_xz);
		
		fftw_free(Txx_phys);
		fftw_free(Tyy_phys);
		fftw_free(Tzz_phys);
		fftw_free(Txy_phys);
		fftw_free(Txx_xz);
		fftw_free(Tyy_xz);
		fftw_free(Tzz_xz);
		fftw_free(Txy_xz);
	}
	
	
	return;
}





extern void dump_reystress_les(
                               double *u_xz,
                               double *v_xz,
                               double *w_xz,
                               double *Txx_xz,
                               double *Tyy_xz,
                               double *Tzz_xz,
                               double *Txy_xz,
                               parameters *par)
{
	int j;
    double *uu, *vv, *ww, *uv;
	double *Txx, *Tyy, *Tzz, *Txy;
	double *ratio;
	
    const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
	
	uu = dvector(0, ny-1);
	vv = dvector(0, ny-1);
	ww = dvector(0, ny-1);
	uv = dvector(0, ny-1);
	Txx = dvector(0, ny-1);
	Tyy = dvector(0, ny-1);
	Tzz = dvector(0, ny-1);
	Txy = dvector(0, ny-1);
	ratio = dvector(0, ny-1);
	
	
	
	// Initialise to zero.
    for (j=0; j<=ny-1; ++j){
		uu[j] = 0.0;
		vv[j] = 0.0;
		ww[j] = 0.0;
		uv[j] = 0.0;
		Txx[j] = 0.0;
		Tyy[j] = 0.0;
		Tzz[j] = 0.0;
		Txy[j] = 0.0;
		ratio[j] = 0.0;
	}
	
	
	profiles_reystress(u_xz, u_xz, uu, par);
	profiles_reystress(v_xz, v_xz, vv, par);
	profiles_reystress(w_xz, w_xz, ww, par);
	profiles_reystress(u_xz, v_xz, uv, par);
	
	
	// DC conponent i.e. average of Tij
	int i, k, index;
	double sign = 1.0;
	i = 0;
	k = 0;
	for (j=0; j<=ny-1; ++j){
		index = 2*k + i*(2*(nz/2+1)) + j*(nx*(2*(nz/2+1)));
        
		if ( (j+0.5)*par->dy < 1.0) sign = - 1.0;
		
		Txx[j] = sqrt(Txx_xz[index]*Txx_xz[index]+Txx_xz[index+1]*Txx_xz[index+1]);
		Tyy[j] = sqrt(Tyy_xz[index]*Tyy_xz[index]+Tyy_xz[index+1]*Tyy_xz[index+1]);
		Tzz[j] = sqrt(Tzz_xz[index]*Tzz_xz[index]+Tzz_xz[index+1]*Tzz_xz[index+1]);
		Txy[j] = sign*sqrt(Txy_xz[index]*Txy_xz[index]+Txy_xz[index+1]*Txy_xz[index+1]);
		
		uu[j] += Txx[j];
		vv[j] += Tyy[j];
		ww[j] += Tzz[j];
		uv[j] += Txy[j];
		
	}
	
	
	double K_res, K_sgs;
	for (j=0; j<=ny-1; ++j){
		K_res = 1.0/2.0*(uu[j]+vv[j]+ww[j]);
		K_sgs = 1.0/2.0*(Txx[j]+Tyy[j]+Tzz[j]);
		ratio[j] = K_sgs/(K_res+K_sgs);
	}
	
	FILE *fp;
	char filename[128];
	double y;
	
	sprintf(filename, "./outputdir/reystress_les_it%d.dat", par->it);
	fp = fopen( filename, "w" );
	fprintf( fp, "variables=j,y,yre,uu,vv,ww,uv,Txx,Tyy,Tzz,Txy,sgs_ratio\n");  	
	for ( j = 0; j <= ny/2; ++j ) {
		y = (j+0.5)*par->dy;
		fprintf( fp, "%d,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e\n",
				j,
				y,
				y*par->re,
				uu[j],
				vv[j],
				ww[j],
				-uv[j],
				Txx[j],
				Tyy[j],
				Tzz[j],
				-Txy[j],
				ratio[j]);
	}
	
	fclose( fp );
	
	
	sprintf(filename, "./outputdir/reystress_les_it_top%d.dat",par->it);
	fp = fopen( filename, "w" );
	fprintf( fp, "variables=j,y,yre,uu,vv,ww,Txx,Tyy,Tzz,Txy,ratio\n");  		
	for ( j = 0; j <= ny/2; ++j ) {
		y = (j+0.5)*par->dy;
		fprintf( fp,"%d,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e\n",
				j,
				y,
				y*par->re,
				uu[ny-1-j],
				vv[ny-1-j],
				ww[ny-1-j],
				uv[ny-1-j],
				Txx[ny-1-j],
				Tyy[ny-1-j],
				Tzz[ny-1-j],
				Txy[ny-1-j],
				ratio[ny-1-j]);
	}
	
	fclose( fp );
	
	free_dvector(uu, 0, ny-1);
	free_dvector(vv, 0, ny-1);
	free_dvector(ww, 0, ny-1);
	free_dvector(uv, 0, ny-1);
	
	free_dvector(Txx, 0, ny-1);
	free_dvector(Tyy, 0, ny-1);
	free_dvector(Tzz, 0, ny-1);
	free_dvector(Txy, 0, ny-1);
	free_dvector(ratio, 0, ny-1);
	
	
	return;
}



extern void statistics(
                       interpolated *inter,
                       fftwplans *planptr,
                       parameters *par )
{
	
	int	i, j, k;
	int	nx, ny, nz;
    
	nx = par->nx;
	ny = par->ny;
	nz = par->nz;
	
	int local_size_wp = par->local_size_wp;
	int	start = par->local_nx_start;
	int end = par->local_nx_start + par->local_nx;
	
	double *scratch;
	double *u_xz;
	double *v_xz;
	double *w_xz;
	double *u_phys;
	double *v_phys;
	double *w_phys;
	
	
	scratch = FFTW_MALLOC(par->local_size_wp);
//if(par->my_rank ==0){
//printf("statistics/ -------------------------get u_physical data defined at center begin\n ");
//}	
	if (par->my_rank == 0){
		u_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		v_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		w_xz = FFTW_MALLOC(nx*ny*2*(nz/2+1));
		u_phys = FFTW_MALLOC(nx*ny*nz);
		v_phys = FFTW_MALLOC(nx*ny*nz);
		w_phys = FFTW_MALLOC(nx*ny*nz);
        
		for (i=0; i<=nx*ny*2*(nz/2+1)-1; ++i){
			u_xz[i] = 0.0;
			v_xz[i] = 0.0;
			w_xz[i] = 0.0;
		}
		
	}else{
		u_xz = NULL;
		v_xz = NULL;
		w_xz = NULL;
		u_phys =NULL;
		v_phys = NULL;
		w_phys = NULL;
	}
//if(par->my_rank ==0){
//printf("statistics/ get u_physical data defined at center begin\n ");
//}	
	// get physical data defined at center
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (par->fd_order == 4){
					scratch[par->index_wp[i][j][k]] \
                    = 9.0/8.0*inter->u_phys_ct1[par->index_ct[i][j][k]] \
                    - 1.0/8.0*inter->u_phys_ct3[par->index_ct[i][j][k]];
				}
				if (par->fd_order == 2){
					scratch[par->index_wp[i][j][k]] \
                    = inter->u_phys_ct1[par->index_ct[i][j][k]];
				}
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
                u_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
 
//if(par->my_rank ==0){
//printf("statistics/ get v_physical data defined at center begin\n ");
//}	   
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				if (par->fd_order == 4){
					scratch[par->index_wp[i][j][k]] \
                    = 9.0/8.0*inter->v_phys_ct1[par->index_ct[i][j][k]] \
                    - 1.0/8.0*inter->v_phys_ct3[par->index_ct[i][j][k]];
                    
				}
				if (par->fd_order == 2){
					scratch[par->index_wp[i][j][k]] \
                    = inter->v_phys_ct1[par->index_ct[i][j][k]];
				}
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
                v_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
//if(par->my_rank ==0){
//printf("statistics/ get w_physical data defined at center begin\n ");
//}	  
	
	for (i=start; i<=end-1; ++i){
		for (j=0; j<=ny-1; ++j){
			for (k=0; k<=nz-1; ++k){
				scratch[par->index_wp[i][j][k]] \
                = inter->w_phys_ct[par->index_ct[i][j][k]];
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    
	MPI_Gatherv( scratch, local_size_wp,             MPI_DOUBLE,
                w_phys,  par->recvcounts_wp, par->displs_wp, MPI_DOUBLE,
                0, MPI_COMM_WORLD );
	
	
	fftw_free(scratch);

//if(par->my_rank ==0){
//printf("statistics/ get w_physical data defined at center begin\n ");
//}    
    
    
	if (par->my_rank == 0){
		
		forward_xz(u_phys, u_xz, planptr, par);
		forward_xz(v_phys, v_xz, planptr, par);
		forward_xz(w_phys, w_xz, planptr, par);
//printf("statistics/ forward_xz succeed\n ");        
		// get reystress
		dump_reystress(u_xz, v_xz, w_xz, par);
		
        
	}

//if(par->my_rank ==0){
//printf("statistics/ dump_reystress succeed\n ");
//}	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (par->my_rank == 0){
		fftw_free(u_phys);
		fftw_free(v_phys);
		fftw_free(w_phys);
		fftw_free(u_xz);
		fftw_free(v_xz);
		fftw_free(w_xz);
	}
	
    
	return;
}





extern void dump_reystress(
                           double *u_xz,
                           double *v_xz,
                           double *w_xz,
                           parameters *par)
{
	int j;
    int nx, nz, ny;
    double *uu, *vv, *ww;
	double *uv, *uw, *vw;
    
    nx = par->nx;
    nz = par->nz;
    ny = par->ny;
    
	uu = dvector(0, ny-1);
	vv = dvector(0, ny-1);
	ww = dvector(0, ny-1);
	
	uv = dvector(0, ny-1);
	uw = dvector(0, ny-1);
	vw = dvector(0, ny-1);
	
	
    // Initialise to zero.
    for (j=0; j<=ny-1; ++j){
		uu[j] = 0.0;
		vv[j] = 0.0;
		ww[j] = 0.0;
		
		uv[j] = 0.0;
		uw[j] = 0.0;
		vw[j] = 0.0;
	}
	
	profiles_reystress(u_xz, u_xz, uu, par);
	profiles_reystress(v_xz, v_xz, vv, par);
	profiles_reystress(w_xz, w_xz, ww, par);
	
	profiles_reystress(u_xz, v_xz, uv, par);
	profiles_reystress(v_xz, w_xz, vw, par);
	profiles_reystress(u_xz, w_xz, uw, par);
	
	
	FILE *fp;
	char filename[128];
	double y;
    
	sprintf(filename, "./outputdir/reystress_it%d.dat", par->it);  //change on 9.27
	fp = fopen( filename, "w" );
	fprintf( fp, "variables=j,y,yre,uu,vv,ww,uv,uw,vw,uu1,vv1,ww1,uv1,uw1,vw1\n");  
    for ( j = 0; j <= ny/2; ++j ) {
		y = (j+0.5)*par->dy;
		fprintf( fp, "%d,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e,%+e\n",
                j,
                y,
                y*par->re,
                uu[j],
                vv[j],
                ww[j],
                uv[j],
                uw[j],
                vw[j],
                uu[ny-1-j],
                vv[ny-1-j],
                ww[ny-1-j],
                uv[ny-1-j],
                uw[ny-1-j],
                vw[ny-1-j]);
	}
	
	fclose( fp );
	
	free_dvector(uu, 0, ny-1);
	free_dvector(vv, 0, ny-1);
	free_dvector(ww, 0, ny-1);
	
	free_dvector(uv, 0, ny-1);
	free_dvector(uw, 0, ny-1);
	free_dvector(vw, 0, ny-1);
	
	return;
}



extern void profiles_reystress(
                               double *u_xz,
                               double *v_xz,
                               double *uv,
                               parameters *par )
{
	int i, j, k;
	int index;
    int nx, nz, ny;
    double norm;
    
    nx = par->nx;
    nz = par->nz;
    ny = par->ny;
    
	for (j=0; j<=ny-1; ++j){
		for (i=0; i<=nx-1; ++i){
			for (k=0; k<=(nz/2+1)-1; ++k){
                
				index = 2*k + i*(2*(nz/2+1)) + j*(nx*(2*(nz/2+1)));
				
                norm = u_xz[index]  *v_xz[index]
                + u_xz[index+1]*v_xz[index+1];
                
                if ( k==0 || k==nz/2 ) {
                    // Minus DC component.
                    if (!( i == 0 && k == 0 )){
                        uv[j] += norm;
					}
				}else{
                    // Contributions from k and -k wavenumbers
                    uv[j] += 2.0*norm;
				}
            }
		}
	}
	
	return;
}



void forward_xz(
                double *u_phys,
                double *u_xz,
                fftwplans *planptr,
                parameters *par)
{
	int i, j, k, index;
	int nx, ny, nz;
    
	nx = par->nx;
	ny = par->ny;
	nz = par->nz;
	
	double *scratch1;
	
	scratch1 = FFTW_MALLOC(nx*ny*2*(nz/2+1));
	
	index = 0;
	for (j=0; j<=ny-1; ++j){
        for (i=0; i<=nx-1; ++i){
            for (k=0; k<=2*(nz/2+1)-1; ++k){
				if (k<=nz-1){
					scratch1[index] = u_phys[par->index_bigwp[i][j][k]];
				}else{
					scratch1[index] = 0.0;
				}
				index = index + 1;
				
			}
		}
	}
	
	fftw_execute_dft_r2c(planptr->p2d_xz_wp, scratch1, (fftw_complex*)u_xz);
	
	for (i=0; i<=nx*ny*2*(nz/2+1)-1; ++i){
		u_xz[i] = u_xz[i]/(double)(nx*nz);
	}
	
	fftw_free(scratch1);
	
	return;
}




extern void courant(
                    double *u_phys,
                    double *v_phys,
                    double *w_phys,
                    parameters *par )
{
    int i;
    double dt, dx, dy, dz;
    double local_courant_no, courant_no;
    double point_courant_no_x, point_courant_no_y, point_courant_no_z;
    
    dt = par->dt;
    dx = par->dx;
    dy = par->dy;
	dz = par->dz;
    
    local_courant_no = 0.0;
    
	for (i=0; i<=par->local_size_u-1; ++i){
		point_courant_no_x = fabs(u_phys[i])*dt/dx;
        
		if ( point_courant_no_x > local_courant_no ){
			local_courant_no = point_courant_no_x;
		}
	}
	
	for (i=0; i<=par->local_size_v-1; ++i){
		point_courant_no_y = fabs(v_phys[i])*dt/dy;
        
		if ( point_courant_no_y > local_courant_no ){
			local_courant_no = point_courant_no_y;
		}
	}
	
	for (i=0; i<=par->local_size_wp-1; ++i){
		point_courant_no_z = fabs(w_phys[i])*dt/dz;
        
		if ( point_courant_no_z > local_courant_no ){
			local_courant_no = point_courant_no_z;
		}
	}
	
    MPI_Reduce( &local_courant_no, &courant_no,
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    
    if ( par->my_rank == 0 ) {
        
		par->cfl = courant_no;
        
		if (par->it%par->choice == 0 ){
			FILE *fp;
			char filename[100];
			
			sprintf( filename, "./outputdir/courant.dat");
            
			fp = fopen( filename, "a" );
            
			fprintf( fp, "%d %+e %+e %+e\n",
                    par->it, par->tott, par->dt, par->cfl );
            
			fclose( fp );
		}
    }
    
    MPI_Bcast( &par->cfl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
}
