/*
 *  correction.c
 * 3D channel flow
 * Periodic in x- and z- direction
 *
 * The last one of the three fractional steps.
 * Calculates the pressure gradients and adds them to the intermediate velocities.
 */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "derivatives__2.h"
#include "derivatives__4.h"
#include "correction.h"


void get_correction(
                    variables *var1,
                    int substep,
                    parameters *par)
{
    
	const double dt = par->dt;
	int i;
	int tag;
	const int local_size_u = par->local_size_u;
	const int local_size_v = par->local_size_v;
	const int local_size_wp = par->local_size_wp;
	
	const double alpha = par->rk3_alpha[substep];
	const double beta = par->rk3_beta[substep];
	
	tag = 300;
	if (par->fd_order == 2){
		get_gradient__2(var1, tag, par);
	}
	if (par->fd_order == 4){
		get_gradient__4(var1, tag, par);
	}
    
	for (i=0; i<=local_size_u-1; ++i){
		var1->u[i] = var1->u[i] - dt*(alpha+beta)*var1->Gpx[i];
	}
	for (i=0; i<=local_size_v-1; ++i){
		var1->v[i] = var1->v[i] - dt*(alpha+beta)*var1->Gpy[i];
	}
	for (i=0; i<=local_size_wp-1; ++i){
		var1->w[i] = var1->w[i] - dt*(alpha+beta)*var1->Gpz[i];
	}
    
	return;
}
