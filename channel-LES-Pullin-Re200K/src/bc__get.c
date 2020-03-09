/*
 * bc__get.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * Implements the log-relation for the wall model.
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "float.h"
#include "definitions.h"
#include "nrutil.h"
#include "bc__get.h"

#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )


/* ******************************************************************* */
// For the initial conditions, bc_dudy is needed.
// bc_u and bc_K are obtained from the initial velocity field.
/* ******************************************************************* */
extern void bc__get_bc_dudy(
							double *bc_u,
							double *bc_K,
							double *bc_dudy,
							char *kinds,
							parameters *par )
{
    int i, k;
	int index;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double nu = 1.0/par->re;
    const double h0 = par->h0;
	double K, u, u_tau, sgn_u;
	
	double bc_sign = 0.0;
	
	if (strcmp(kinds, "bot")==0){
		bc_sign = 1.0;
	}else if (strcmp(kinds, "top")==0){
		bc_sign = - 1.0;
	}else{
		printf("error bc__get_bc_dudy");
	}
	
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
            
			sgn_u = bc_u[index] > 0.0 ? 1.0 : -1.0;
			u = fabs( bc_u[index] );
			K = bc_K[index];
			K = 0.41;
			u_tau = bc__u_tau_from_u( u, 0.0, K, 11.0, nu, h0);
			bc_dudy[index] = bc_sign*sgn_u*u_tau*u_tau/nu;
		}
	}
	
}


/* ******************************************************************* */
// For the initial conditions, bc_dudy is needed.
// Get u_tau from the initial velocty field.
// bc_u and bc_K are obtained from the initial velocity field.
/* ******************************************************************* */
extern double bc__u_tau_from_u(
							   double u,
							   double gamma,
							   double K,
							   double z0p,
							   double nu,
							   double z )
{
    double u_tau, u_tau_old, z_plus, f, dfdu_tau, kappa, dkappadu_tau;
    const double omega = 0.5; // Under-relaxation
	
    // The linear part. ( u/u_tau = z u_tau/nu ).
    u_tau = sqrt( u*nu/z );
    z_plus = z*u_tau/nu;
	
    // The log part.
    if ( z_plus > z0p ) {
        // Newton iterations.
        do {
            kappa = K;
            dkappadu_tau = 0.0;
            f = u/u_tau - log( z_plus/z0p )/kappa - z0p;
            dfdu_tau = ( - u/u_tau - 1.0/kappa )/u_tau
			+ log( z_plus/z0p )*dkappadu_tau/( kappa*kappa );
            u_tau_old = u_tau;
            u_tau = u_tau_old - omega*f/dfdu_tau;
            if ( u_tau < 0.0 ) u_tau = FLT_EPSILON;
            z_plus = z*u_tau/nu;
        }
        while ( fabs( 1.0 - u_tau_old/u_tau ) > FLT_EPSILON );
    }
	
    return u_tau;
}

/* ******************************************************************* */
// bc_u is calculated.
/* ******************************************************************* */
extern void bc__get_bc_vel(
						   double *bc_dudy,
						   double *bc_K,
						   double *bc_kappa,
						   double *bc_u,
						   double *bc_v,
						   double *bc_w,
						   char *kinds,
						   parameters *par )
{
	int i, k;
	int index;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
	const double nu = 1.0/par->re;
    const double h0 = par->h0;
	
	double bc_sign = 0.0;
    double dudy, u_tau, sgn_u, kappa;
	
	if (strcmp(kinds, "bot")==0){
		bc_sign = 1.0;
	}else if (strcmp(kinds, "top")==0){
		bc_sign = - 1.0;
	}else{
		printf("error bc__get_bc_vel");
	}
	
	for (i=start; i<=end-1; ++i) {
		for (k=0; k<=nz-1; ++k) {
			index = par->index_tb[i][k];
			
			dudy = bc_dudy[index];
			u_tau = sqrt( nu*fabs( dudy ) );
            
			kappa = bc_kappa[index];
            
			
			sgn_u = dudy/bc_sign > 0.0 ? 1.0 : -1.0;
			
			bc_u[index]
            = sgn_u*bc__u_from_u_tau( u_tau, 0.0, 0.0, kappa, 11.0, nu, h0 , par);
			bc_v[index] = 0.0;
			bc_w[index] = 0.0;
		}
	}
	
	return;
}

/* ******************************************************************* */
/* Log-law for smooth-walled boundary layer and Colebrook's formula as a
 * roughness correction.
 * The slip velocity is calculated at each wall point and at every time step.
 */
/* ******************************************************************* */
extern double bc__u_from_u_tau(
							   double u_tau,
							   double gamma,
							   double K,
							   double kappa,
							   double z0p,
							   double nu,
							   double z ,
                               parameters *par )

{
	double z_plus = z*u_tau/nu;
	return u_tau*( log( z_plus/z0p/(1.0+0.26*par->eps*u_tau/nu) )/kappa + z0p );
}

