/*
 *  les.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * Tij calculations using the stretched vortex model.
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "definitions.h"
#include "nrutil.h"
#include "convective.h"
#include "les.h"
#include "spiral.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

extern void les_channel(
						variables *var,
						interpolated *inter,
						parameters *par,
						fftwplans *planptr)
{
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
    const int local_nx_near = 1;
	const int ny_near = 1;
	const int nz_near = 1;
    const int size_near
	= ( 2*local_nx_near + 1 )*( 2*ny_near + 1 )*( 2*nz_near + 1 ); // = 27
    int i, j, k;
	int index;
    int i_near, j_near, k_near;
    int p, q, r;
    double dudx[3][3], e[3];
    double *u_near, *v_near, *w_near, *x_near, *y_near, *z_near;
    int *yes_near; // Indicator
    int ix_near, ix_this = 100; // ix_this has to be non zero
    double txx, tyy, tzz, txy, tyz, tzx, K;
    double lv;
    const double nu = 1.0/par->re;
    const double delta = pow( ( par->dx )
                             *( par->dy )
                             *( par->dz ), 1.0/3.0 );
	
    u_near = dvector( 0, size_near - 1 ); //freed
    v_near = dvector( 0, size_near - 1 ); //freed
    w_near = dvector( 0, size_near - 1 ); //freed
    x_near = dvector( 0, size_near - 1 ); //freed
    y_near = dvector( 0, size_near - 1 ); //freed
    z_near = dvector( 0, size_near - 1 ); //freed
    yes_near = ivector( 0, size_near - 1 ); //freed
    

    for ( i=start; i<=end-1; ++i ) {
        for ( j=0; j<=ny; ++j ) {
            for ( k=0; k<=nz-1; ++k) {
				
                // Obtain velocities of neighbouring points.
                ix_near = 0;
                for ( p = - local_nx_near; p <= local_nx_near; p++ ) {
                    for ( q = - ny_near; q <= ny_near; q++ ) {
                        for ( r = - nz_near; r <= nz_near; r++ ) {
							
                            i_near = i + p;
                            j_near = j + q;
				 k_near = k + r;
							
                            if ( j_near < 0 || ny < j_near ) j_near = j - q;
							
							x_near[ix_near] = par->dx*i_near;
                            			y_near[ix_near] = par->dy*j_near;
                            			z_near[ix_near] = par->dz*k_near;
							
							if ( k_near < 0) k_near = nz-1;
							if ( nz - 1 < k_near) k_near = 0;
							
							index = par->index_les[i_near][j_near][k_near];
							
							u_near[ix_near] = inter->u_phys_les[index];
							v_near[ix_near] = inter->v_phys_les[index];
							w_near[ix_near] = inter->w_phys_les[index];
							
                            if ( p == 0 && q == 0 && r == 0 ) {
							ix_this = ix_near;
                            }
							yes_near[ix_near] = 1;
							
                            if ( ix_near == ix_this ) {
								yes_near[ix_near] = 0;
							}
							
                            ix_near++;
                  }
			}
				}
				
//   printf("pqr succeed !i=%d, j=%d,k=%d\n",i,j,k);             
				index = par->index_cn[i][j][k];
				dudx[0][0] = var->dudx[index];
                dudx[0][1] = var->dudy[index];
                dudx[0][2] = var->dudz[index];
                dudx[1][0] = var->dvdx[index];
                dudx[1][1] = var->dvdy[index];
                dudx[1][2] = var->dvdz[index];
                dudx[2][0] = var->dwdx[index];
                dudx[2][1] = var->dwdy[index];
                dudx[2][2] = var->dwdz[index];
				
                
                e[0] = 0.0;
                e[1] = 0.0;
                e[2] = 0.0;
                
				
				spiral_sgs_stress(u_near, v_near, w_near,
								  x_near, y_near, z_near, yes_near, size_near, ix_this,
								  dudx, e, nu, delta,
								  &txx, &tyy, &tzz,
								  &txy, &tyz, &tzx, &K, &lv);
//   printf("spiral_sgs_stress succeed !i=%d, j=%d,k=%d\n",i,j,k);      				
				index = par->index_cn[i][j][k];
                var->Txx[index] = txx;
                var->Tyy[index] = tyy;
                var->Tzz[index] = tzz;
                var->Txy[index] = txy;
                var->Tyz[index] = tyz;
                var->Tzx[index] = tzx;
                var->K[index] = K;
                
                
                index = par->index_tb[i][k];
				if (j == 0) {
					var->bc_txy_b[index] = txy;
				}
				if (j == ny) {
					var->bc_txy_t[index] = txy;
				}
				
				
                if ( (j == 1 || j == ny-1) && par->bc == 1) { // next to the wall
					
					
                    
                    
					spiral_sgs_stress(u_near, v_near, w_near,
									  x_near, y_near, z_near, yes_near, size_near, ix_this,
									  dudx, e, 0.0, delta,
									  &txx, &tyy, &tzz,
									  &txy, &tyz, &tzx, &K, &lv);
					
                    
					// Streamwise vortices
					e[0] = 1.0;
					e[1] = 0.0;
					e[2] = 0.0;
					
					
					spiral_sgs_stress(u_near, v_near, w_near,
									  x_near, y_near, z_near, yes_near, size_near, ix_this,
									  dudx, e, 0.0, delta,
									  &txx, &txx, &txx,
									  &txx, &txx, &txx, &K, &lv);
                    
					index = par->index_tb[i][k];
                    if (j == 1) {
						var->bc_kappa_b[index] = 0.45/2.0*sqrt(K)/sqrt(fabs(txy));
						var->bc_K_b[index] = K;
					}
					if (j == ny-1) {
						var->bc_kappa_t[index] = 0.45/2.0*sqrt(K)/sqrt(fabs(txy));
						var->bc_K_t[index] = K;
					}
					
					
                }
			}
		}
	}
if(par->my_rank == 0){
// ("spiral_sgs_stress succeed---------------------------------------------------------------\n");
}		
    free_dvector( u_near, 0, size_near - 1 );
    free_dvector( v_near, 0, size_near - 1 );
    free_dvector( w_near, 0, size_near - 1 );
    free_dvector( x_near, 0, size_near - 1 );
    free_dvector( y_near, 0, size_near - 1 );
    free_dvector( z_near, 0, size_near - 1 );
    free_ivector( yes_near, 0, size_near - 1 );
	
	
	
	fftw_execute_r2r(planptr->p1d_z_cn, var->Txx, var->Txx);
	fftw_execute_r2r(planptr->p1d_z_cn, var->Tyy, var->Tyy);
	fftw_execute_r2r(planptr->p1d_z_cn, var->Tzz, var->Tzz);
	fftw_execute_r2r(planptr->p1d_z_cn, var->Txy, var->Txy);
	fftw_execute_r2r(planptr->p1d_z_cn, var->Tyz, var->Tyz);
	fftw_execute_r2r(planptr->p1d_z_cn, var->Tzx, var->Tzx);
	for (i=0; i<=par->local_size_cn-1; ++i){
		var->Txx[i] = var->Txx[i]/(double)(nz);
		var->Tyy[i] = var->Tyy[i]/(double)(nz);
		var->Tzz[i] = var->Tzz[i]/(double)(nz);
		var->Txy[i] = var->Txy[i]/(double)(nz);
		var->Tyz[i] = var->Tyz[i]/(double)(nz);
		var->Tzx[i] = var->Tzx[i]/(double)(nz);
	}
    
    
    if (par->fd_order == 2){
		trancate_exp(var->Txx, "cn2", par);
		trancate_exp(var->Tyy, "cn2", par);
		trancate_exp(var->Tzz, "cn2", par);
		trancate_exp(var->Txy, "cn2", par);
		trancate_exp(var->Tyz, "cn2", par);
		trancate_exp(var->Tzx, "cn2", par);
	}
	if (par->fd_order == 4){
		trancate_exp(var->Txx, "cn4", par);
		trancate_exp(var->Tyy, "cn4", par);
		trancate_exp(var->Tzz, "cn4", par);
		trancate_exp(var->Txy, "cn4", par);
		trancate_exp(var->Tyz, "cn4", par);
		trancate_exp(var->Tzx, "cn4", par);
	}
 //if(par->my_rank == 0){
//printf("trancate_exp succeed\n");
//}   
    return;
}



