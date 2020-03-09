/*
 *  fft.c
 * 3D channel flow
 * Periodic in x- and z- direction
 * 
 * fft related functions are defined here.
 * void fft_make_plans(variables *var, parameters *par, fftwplans *planptr)
 * void fft_make_plans__transpose(parameters *par, fftwplans *planptr)
 * extern void fftw_forward_poisson_ver2(double *phys, double *wave, fftwplans *planptr, parameters *par)
 * extern void fftw_inverse_poisson_ver2(double *wave, double *phys, fftwplans *planptr, parameters *par)
 * void fft_destroy_plans__transpose(fftwplans *planptr, parameters *par)
 * void fft_destroy_plans(fftwplans *planptr, parameters *par)
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "fftw3.h"
#include "fftw3-mpi.h"
#include "definitions.h"
#include "nrutil.h"
#include "parameters.h"
#include "fft.h"
#include "mpi.h"





void fft_make_plans(variables *var, parameters *par, fftwplans *planptr)
{
	int i;
	const int local_nx = par->local_nx;
	const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
	
	const int local_size_ued = par->local_size_ued;
	const int local_size_ved = par->local_size_ved;
	const int local_size_ct = par->local_size_ct;
	const int local_size_cn = par->local_size_cn;
	
	fftw_r2r_kind	kind_f; kind_f = FFTW_R2HC;
	fftw_r2r_kind	kind_i; kind_i = FFTW_HC2R;
	
	double *sample_ued, *sample_ved, *sample_ct, *sample_cn;
	
	sample_ued = FFTW_MALLOC(local_size_ued);
	sample_ved = FFTW_MALLOC(local_size_ved);
	sample_ct = FFTW_MALLOC(local_size_ct);
	sample_cn = FFTW_MALLOC(local_size_cn);
	
	for (i=0; i<=local_size_ued-1; ++i){
		sample_ued[i] = 0.0;
	}
	for (i=0; i<=local_size_ved-1; ++i){
		sample_ved[i] = 0.0;
	}
	for (i=0; i<=local_size_ct-1; ++i){
		sample_ct[i] = 0.0;
	}
	for (i=0; i<=local_size_cn-1; ++i){
		sample_cn[i] = 0.0;
	}
    
	// forward transform plan for velocity u
	planptr->p1d_z_u = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ local_nx*ny, \
                                          var->u, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                          var->u, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                          /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	
	
	// forward transform plan for velocity v
	planptr->p1d_z_v = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(ny-1)*local_nx, \
                                          var->v, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                          var->v, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                          /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	
	// forward transform plan for velocity w or p
	planptr->p1d_z_wp = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ local_nx*ny, \
                                           var->w, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                           var->w, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                           /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	
    
	// forward transform plan for top and bottom boundary
	planptr->p1d_z_bc_tb = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ local_nx+1, \
                                              var->bc_u_t, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                              var->bc_u_t, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                              /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
    
	// forward transform plan for left and right boundary
	if (par->my_rank == 0){
		planptr->p1d_z_bc_uw_lr = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ ny+2, \
                                                     var->bc_u_l, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                     var->bc_u_l, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                     /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
		planptr->p1d_z_bc_v_lr = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ ny+1, \
                                                    var->bc_v_l, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                    var->bc_v_l, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                    /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	}
	
	// forward transform plan for left and right boundary
	if (par->my_rank == par->nproc - 1){
		planptr->p1d_z_bc_uw_lr = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ ny+2, \
                                                     var->bc_u_r, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                     var->bc_u_r, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                     /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
		planptr->p1d_z_bc_v_lr = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ ny+1, \
                                                    var->bc_v_r, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                    var->bc_v_r, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                    /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	}
    
	
	
	if (par->fd_order == 2){
		// forward transform plan for velocity ued
		planptr->p1d_z_ued = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+2)*(ny+2), \
                                                sample_ued, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                sample_ued, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ved
		planptr->p1d_z_ved = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+2)*(ny+1), \
                                                sample_ved, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                sample_ved, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ct
		planptr->p1d_z_ct = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+2)*(ny+2), \
                                               sample_ct, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                               sample_ct, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                               /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity cn
		planptr->p1d_z_cn = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+1)*(ny+1), \
                                               sample_cn, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                               sample_cn, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                               /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	}
	
	
	
	if (par->fd_order == 4){
		// forward transform plan for velocity ued
		planptr->p1d_z_ued = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+6)*(ny+4), \
                                                sample_ued, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                sample_ued, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ved
		planptr->p1d_z_ved = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+6)*(ny+5), \
                                                sample_ved, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                sample_ved, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ct
		planptr->p1d_z_ct = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+6)*(ny+4), \
                                               sample_ct, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                               sample_ct, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                               /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity cn
		planptr->p1d_z_cn = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+3)*(ny+3), \
                                               sample_cn, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                               sample_cn, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                               /*kind*/&kind_f, /*flag*/FFTW_PATIENT);
	}
    
	
	
	
	
	
	
	
	
	
	
	
	
	
	// inverse transform plan for velocity u
	planptr->p1d_invz_u = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ local_nx*ny, \
                                             var->u, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                             var->u, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                             /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
    
	// inverse transform plan for velocity v
	planptr->p1d_invz_v = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(ny-1)*local_nx, \
                                             var->v, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                             var->v, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                             /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
    
	// inverse transform plan for velocity w or p
	planptr->p1d_invz_wp = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ local_nx*ny, \
                                              var->w, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                              var->w, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                              /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
    
	// inverse transform plan for top and bottom boundary
	planptr->p1d_invz_bc_tb = fftw_plan_many_r2r(/*rank*/1,  &nz, /*howmany*/ local_nx+1, \
                                                 var->bc_u_t, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                 var->bc_u_t, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                 /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
    
	// inverse transform plan for left and right boundary
	if ( (par->my_rank == 0) ){
		planptr->p1d_invz_bc_uw_lr = fftw_plan_many_r2r(/*rank*/1,  &nz, /*howmany*/ ny+2, \
                                                        var->bc_u_l, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                        var->bc_u_l, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                        /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
		planptr->p1d_invz_bc_v_lr = fftw_plan_many_r2r(/*rank*/1,  &nz, /*howmany*/ ny+1, \
                                                       var->bc_u_l, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                       var->bc_u_l, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                       /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
	}
	if ( (par->my_rank == par->nproc - 1) ){
		planptr->p1d_invz_bc_uw_lr = fftw_plan_many_r2r(/*rank*/1,  &nz, /*howmany*/ ny+2, \
                                                        var->bc_u_r, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                        var->bc_u_r, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                        /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
		planptr->p1d_invz_bc_v_lr = fftw_plan_many_r2r(/*rank*/1,  &nz, /*howmany*/ ny+1, \
                                                       var->bc_u_r, /*inembed*/ &nz, /*istride*/ 1, /*idist*/ nz, \
                                                       var->bc_u_r, /*onembed*/ &nz, /*ostride*/ 1, /*odist*/ nz, \
                                                       /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
	}
    
    
    
	
	if (par->fd_order == 2){
		// forward transform plan for velocity ued
		planptr->p1d_invz_ued = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+2)*(ny+2), \
                                                   sample_ued, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                   sample_ued, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                   /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ved
		planptr->p1d_invz_ved = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+2)*(ny+1), \
                                                   sample_ved, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                   sample_ved, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                   /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ct
		planptr->p1d_invz_ct = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+2)*(ny+2), \
                                                  sample_ct, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                  sample_ct, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                  /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity cn
		planptr->p1d_invz_cn = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+1)*(ny+1), \
                                                  sample_cn, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                  sample_cn, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                  /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
	}
	
	if (par->fd_order == 4){
		// forward transform plan for velocity ued
		planptr->p1d_invz_ued = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+6)*(ny+4), \
                                                   sample_ued, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                   sample_ued, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                   /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ved
		planptr->p1d_invz_ved = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/(local_nx+6)*(ny+5), \
                                                   sample_ved, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                   sample_ved, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                   /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity ct
		planptr->p1d_invz_ct = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+6)*(ny+4), \
                                                  sample_ct, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                  sample_ct, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                  /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
        
		// forward transform plan for velocity cn
		planptr->p1d_invz_cn = fftw_plan_many_r2r(/*rank*/1, &nz, /*howmany*/ (local_nx+3)*(ny+3), \
                                                  sample_cn, /*inembed*/&nz, /*istride*/ 1, /*idist*/ nz, \
                                                  sample_cn, /*onembed*/&nz, /*ostride*/ 1, /*odist*/ nz, \
                                                  /*kind*/&kind_i, /*flag*/FFTW_PATIENT);
	}
	
	
	
	
	
	fftw_free(sample_ued);
	fftw_free(sample_ved);
	fftw_free(sample_ct);
	fftw_free(sample_cn);
	
	
	
	double *in;
	double *out;
	int xz_2d[2], inembed[2], onembed[2];
	
	if (par->my_rank == 0){
		in = FFTW_MALLOC(ny*nx*2*(nz/2+1));
		out = FFTW_MALLOC(ny*nx*2*(nz/2+1));
        
		xz_2d[0] = nx;
		xz_2d[1] = nz;
		inembed[0] = nx;
		inembed[1] = 2*(nz/2+1);
		onembed[0] = nx;
		onembed[1] = (nz/2+1);
        
		planptr->p2d_xz_wp = fftw_plan_many_dft_r2c(/*rank*/2, xz_2d, /*howmany*/ ny,
                                                    /*double in*/ in, /*inembed*/inembed, /*istride*/1 , /*idist*/nx*2*(nz/2+1),
                                                    /*complex out*/ (fftw_complex*)out, /*onembed*/onembed, /*ostride*/1 , /*odist*/nx*(nz/2+1) ,
                                                    /*flag*/FFTW_PATIENT);
		
		xz_2d[0] = nx;
		xz_2d[1] = nz;
		inembed[0] = nx;
		inembed[1] = (nz/2+1);
		onembed[0] = nx;
		onembed[1] = 2*(nz/2+1);
        
		planptr->p2d_invxz_wp = fftw_plan_many_dft_c2r(/*rank*/2, xz_2d, /*howmany*/ ny,
                                                       /*complex in*/ (fftw_complex*)in, /*inembed*/inembed, /*istride*/1 , /*idist*/nx*(nz/2+1),
                                                       /*double out*/ out, /*onembed*/onembed, /*ostride*/1 , /*odist*/nx*2*(nz/2+1),
                                                       /*flag*/FFTW_PATIENT);
        
        
		fftw_free(in);
		fftw_free(out);
        
	}
	
	
	return;
}



void fft_make_plans__transpose(parameters *par, fftwplans *planptr)
{
	
	int i;
	const int nx = par->nx;
	const int ny = par->ny;
	const int nz = par->nz;
	const int local_nx = par->local_nx;
	const int local_nz = par->nz/par->nproc;
	
	double *in, *out;
	
    
    {
        in = FFTW_MALLOC(ny*local_nx*nz);
        out = FFTW_MALLOC(ny*nx*local_nz);
        
        for (i=0; i<=ny*local_nx*nz-1; ++i){
            in[i] = 0.0;
        }
        for (i=0; i<=ny*nx*local_nz-1; ++i){
            out[i] = 0.0;
        }
        
        const ptrdiff_t N0 = nx, N1 = nz, howmany = ny;
        const ptrdiff_t block0 = local_nx, block1 = local_nz;
        
        
        planptr->transpose_yz2xy = fftw_mpi_plan_many_transpose(N0, N1, howmany, block0, block1,
                                                                in, out, MPI_COMM_WORLD, FFTW_PATIENT);
        
        planptr->transpose_xy2yz = fftw_mpi_plan_many_transpose(N1, N0, howmany, block1, block0,
                                                                out, in, MPI_COMM_WORLD, FFTW_PATIENT);
        
        
        fftw_free(in);
        fftw_free(out);
    }
    
    
    {
        fftw_r2r_kind	kind_f; kind_f = FFTW_R2HC;
        fftw_r2r_kind	kind_i; kind_i = FFTW_HC2R;
        
        in = FFTW_MALLOC(ny*nx*local_nz);
        out = FFTW_MALLOC(ny*nx*local_nz);
        
        for (i=0; i<=ny*nx*local_nz-1; ++i){
            out[i] = 0.0;
            in[i] = 0.0;
        }
        
        // forward transform plan in x for velocity w or p after transpose
        planptr->p1d_x_wp = fftw_plan_many_r2r(1, &nx, local_nz*ny,
                                               in, &nx, 1, nx,
                                               out, &nx, 1, nx,
                                               &kind_f, FFTW_PATIENT);
        
        planptr->p1d_invx_wp = fftw_plan_many_r2r(1, &nx, local_nz*ny,
                                                  out, &nx, 1, nx,
                                                  in, &nx, 1, nx,
                                                  &kind_i, FFTW_PATIENT);
        
        fftw_free(in);
        fftw_free(out);
    }
    
    {
        fftw_r2r_kind	kind_f = FFTW_REDFT10;
        fftw_r2r_kind	kind_i = FFTW_REDFT01;
        
        // transform plan for velocity wp
        in = FFTW_MALLOC(par->local_size_wp);
        out = FFTW_MALLOC(par->local_size_wp);
        for (i=0; i<=par->local_size_wp-1; ++i){
            in[i] = 0.0;
            out[i] = 0.0;
        }
        
        
        // Cosine transform in y used before data transpose
        // stored in [i][k][j]
        const int howmany = local_nx*nz;
        planptr->p1d_poisson_ver2 = fftw_plan_many_r2r(1, &ny, howmany,
                                                       in, &ny, 1, ny,
                                                       out, &ny, 1, ny,
                                                       &kind_f, FFTW_PATIENT);
        
        planptr->p1d_invpoisson_ver2 = fftw_plan_many_r2r(1, &ny, howmany,
                                                          out, &ny, 1, ny,
                                                          in, &ny, 1, ny,
                                                          &kind_i, FFTW_PATIENT);
        
        
        
        fftw_free(in);
        fftw_free(out);
    }
    
	return;
}


extern void fftw_forward_poisson_ver2(double *phys, double *wave, fftwplans *planptr, parameters *par)
{
	int i;
    const double inv2ny = 1.0/(double)(2*par->ny);
    
	double *copy;
	copy = FFTW_MALLOC(par->local_size_wp);
	for (i=0; i<=par->local_size_wp-1; ++i){
		copy[i] = phys[i];
	}
    
    
	fftw_execute_r2r(planptr->p1d_poisson_ver2, copy, wave);
    
	for (i=0; i<=par->local_size_wp-1; ++i){
		wave[i] = wave[i]*inv2ny;
	}
	fftw_free(copy);
	
	return;
}


extern void fftw_inverse_poisson_ver2(double *wave, double *phys, fftwplans *planptr, parameters *par)
{
	
	int i;
	double *copy;
	copy = FFTW_MALLOC(par->local_size_wp);
	for (i=0; i<=par->local_size_wp-1; ++i){
		copy[i] = wave[i];
	}
	fftw_execute_r2r(planptr->p1d_invpoisson_ver2, copy, phys);
	fftw_free(copy);
	
	
	return;
}


void fft_destroy_plans__transpose(fftwplans *planptr, parameters *par)
{
	
	fftw_destroy_plan(planptr->transpose_yz2xy);
	fftw_destroy_plan(planptr->transpose_xy2yz);
    fftw_destroy_plan(planptr->p1d_x_wp);
    fftw_destroy_plan(planptr->p1d_invx_wp);
    fftw_destroy_plan(planptr->p1d_poisson_ver2);
    fftw_destroy_plan(planptr->p1d_invpoisson_ver2);
	
	return;
}


void fft_destroy_plans(fftwplans *planptr, parameters *par)
{
	//destroying plans
	fftw_destroy_plan(planptr->p1d_z_u);
	fftw_destroy_plan(planptr->p1d_z_v);
	fftw_destroy_plan(planptr->p1d_z_wp);
	fftw_destroy_plan(planptr->p1d_z_bc_tb);
	if ( (par->my_rank == 0) || (par->my_rank == par->nproc - 1)  ){
		fftw_destroy_plan(planptr->p1d_z_bc_uw_lr);
		fftw_destroy_plan(planptr->p1d_z_bc_v_lr);
	}
    
	//destroying plans
	fftw_destroy_plan(planptr->p1d_z_ued);
	fftw_destroy_plan(planptr->p1d_z_ved);
	fftw_destroy_plan(planptr->p1d_z_ct);
	fftw_destroy_plan(planptr->p1d_z_cn);
	
	//destroying plans
	fftw_destroy_plan(planptr->p1d_invz_u);
	fftw_destroy_plan(planptr->p1d_invz_v);
	fftw_destroy_plan(planptr->p1d_invz_wp);
	fftw_destroy_plan(planptr->p1d_invz_bc_tb);
	if ( (par->my_rank == 0) || (par->my_rank == par->nproc - 1) ){
		fftw_destroy_plan(planptr->p1d_invz_bc_uw_lr);
		fftw_destroy_plan(planptr->p1d_invz_bc_v_lr);
	}
	
	//destroying plans
	fftw_destroy_plan(planptr->p1d_invz_ued);
	fftw_destroy_plan(planptr->p1d_invz_ved);
	fftw_destroy_plan(planptr->p1d_invz_ct);
	fftw_destroy_plan(planptr->p1d_invz_cn);
	
	return;
}
