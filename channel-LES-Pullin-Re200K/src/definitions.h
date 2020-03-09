#include <fftw3.h>
#define PI 3.141592653589793
#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )



#ifdef F77_WITH_NO_UNDERSCORE
#define   numroc_      numroc
#define   descinit_    descinit
#define   pdlamch_     pdlamch
#define   pdlange_     pdlange
#define   pdlacpy_     pdlacpy
#define   pdgesv_      pdgesv
#define   pdgemm_      pdgemm
#define   indxg2p_     indxg2p
#endif


extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);
extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, \
						int *icsrc, int *ictxt, int *lld, int *info);
extern double pdlamch_( int *ictxt , char *cmach);
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);

extern void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, \
						double *b, int *ib, int *jb, int *descb);
extern void pddtsv_( int *n, int *nrhs, double *DL, double *D, double *DU, int *ja, int *desca, double *B, int *ib, int *descb, double *work, int *Lwork, int *info);
extern int  indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *D, int *ldb, int *piv, double *scratch, int *N, int *info);
extern void pdgesv_(int *n, int *nrhs, double *DG, int *ia, int *ja, int *descA, int *ippiv, double *x, int *ib, int *jb, int *descB, int *info );





/* ************************************************************************ */
/* Structure for the parameters.                                            */
/* ************************************************************************ */

typedef struct parameters {
        double eps;
	
	int		icont;          //0 random field, 1 read field from file
	int		it;				// iteration
	int		nt;				// Total number of time steps
	int		nx;				// # of grids in x-direction
	int		ny;				// # of grids in y-direction
	int		nz;				// # of grids in z-direction
	int		mz;				// # of grids in z-direction for dealiasing
	int		nq_x;			// order of velocity vector in x-direction
	int		nq_y;			// order of velocity vector in y-direction
	int		nq_z;			// order of velocity vector in w-direction
	
	int local_nx;
	int local_nx_start;
	int local_size_u;
	int local_size_v;
	int local_size_wp;
	
	int local_size_ued;
	int local_size_ved;
	int local_size_ct;
	int local_size_cn;
	
	int local_size_les;
	int local_size_tb;
	
	
	int my_rank;
	int nproc;
	
	int		np;				// order of pressure vector
	
	double	Lx;				// size of box in x-direction
	double	Ly;				// size of box in y-direction
	double	Lz;				// size of box in x-direction

	double	re;				// Reynolds number
	double	re_aim;			// aiming Reynolds number
	double	re_start;		// starting Reynolds number
	double	U;				// mean velocity
	double	dt;				// time step size
	double	cfl;			// curant number
	double	u_tau;			// friction velocity
	double  Du_max;			// max divergence
	
	
	double	tott;			// total time
	int		file_dump;		// frecuency of data output
	int		choice;			//
	int		stat;			//
	int		const_mass;		// constant mass flux 1; constant pressure gradient 0;
	double	mass_flux;		// mass flux
	double	mass_flux_aim;	// aiming mass flux
	double	mass_flux_start;// aiming mass flux
	int		dealias;		// 2/3 rule on 1; 2/3 rule off 0;
	int		non_linear;		// divergence form 1; skew symmetric form 0;
	int		fd_order;		// finite difference scheme order (2 or 4);
	int		width_share;		// finite difference scheme order (2 or 4);
	int		bc_data;

	double	dx;				// grid spacing
	double	dy;
	double	dz;
	
	double	*kz;			// vector of wavenumbers in z direction
	
	int	***index_u;			// index to locate grid position in vector
	int	***index_v;
	int ***index_wp;
	int ***index_bigwp;
	
	
	int ***index_ct;
	int ***index_cn;
	int ***index_ued;
	int ***index_ved;
	
	int ***index_les;

	
	int **index_uw_lr;
	int **index_v_lr;
	int **index_tb;
	
	double	rk3_alpha[3];	// Coefficients for integration scheme
    double	rk3_beta[3];
    double	rk3_gamma[3];
    double	rk3_zeta[3];
	
	double ***inv_DG;
	
	double ***inv_DGx0;
	double ***inv_DGx1;
	double ***inv_DGx2;
	
	double ***inv_DGB_t;
	double ***inv_DGB_b;
	
	int *recvcounts_u;
	int *displs_u;
	int *recvcounts_v;
	int *displs_v;
	int *recvcounts_wp;
	int *displs_wp;
	int *recvcounts_tb;
	int *displs_tb;
    int *recvcounts_mean_2d;
	int *displs_mean_2d;
	
	int *recvcounts_mean_tb;
	int *displs_mean_tb;
	
	char outputdir[100];
	
	int test;
	int les; // SGS model 1: on, 0: off
	int bc; // boundary model 1: on, 0: off
	double manufactured_coeff;
	
	double h0;
	
	double dPdx;
	double mass_flux_curr;
	
	double mean_kappa;
	double mean_dudy;
	double expected_dpdx_t;
	double expected_dpdx_b;
	double mean_bcK_t;
	double mean_bctxy_t;
	double mean_bcK_b;
	double mean_bctxy_b;
	int j_log;
	
	int ode_choice;
	
	double A1;
	double A2;
	double B1;
	double B2;
	double C1;
	double C2;
	double D;
	double E1;
	double E2;
	
} parameters;


/* ************************************************************************ */
/* Structure for the variables                                              */
/* ************************************************************************ */
typedef struct variables {

    double	*u; // velocity vector in x-direction fft-ed in z-direction
    double	*v; // velicity vecotr in y-direction fft-ed in z-direction
	double	*w; // velicity vecotr in x-direction fft-ed in z-direction
	
	double *u_phys; // velicity vecotr in x-direction
    double *v_phys; // velicity vecotr in y-direction
    double *w_phys; // velicity vecotr in z-direction
	
	
	double *dudx; // velocity gradients IN PHYSICAL SPACE
    double *dudy;
    double *dudz;
    double *dvdx;
    double *dvdy;
    double *dvdz;
    double *dwdx;
    double *dwdy;
    double *dwdz;
	
    double *Txx; // subgrid stress defined at cn (par->les == 1)
    double *Tyy; // or ct (par->les == 2)
    double *Tzz;
    double *Txy;
    double *Tyz;
    double *Tzx;
	double *K;
	
	
	double	*Nu; // non-linear term in x momentum 
	double	*Nv; // non-linear term in y momentum 
	double	*Nw; // non-linear term in z momentum 
	
	double	*Lu; // linear term in x momentum 
	double	*Lv; // linear term in y momentum 
	double	*Lw; // linear term in z momentum 
	double	*bc_Lu; // boundary value for linear term in x momentum
	double	*bc_Lv; // boundary value for linear term in x momentum
	double	*bc_Lw; // boundary value for linear term in x momentum
	
	double	*Lu_y; // linear term in x momentum 
	double	*Lv_y; // linear term in y momentum 
	double	*Lw_y; // linear term in z momentum 
	double	*bc_Lu_y; // boundary value for linear term in x momentum
	double	*bc_Lv_y; // boundary value for linear term in x momentum
	double	*bc_Lw_y; // boundary value for linear term in x momentum
	
	double	*Lu_xz; // linear term in x momentum 
	double	*Lv_xz; // linear term in y momentum 
	double	*Lw_xz; // linear term in z momentum 
	double	*bc_Lu_xz; // boundary value for linear term in x momentum
	double	*bc_Lv_xz; // boundary value for linear term in x momentum
	double	*bc_Lw_xz; // boundary value for linear term in x momentum
	
	
	double	*rhs_u; // right hand side in x momentum
	double	*rhs_v; // right hand side in y momentum
	double	*rhs_w; // right hand side in w momentum
	
	double	*bc_Du; // boundary value for divergence i.e. continuity equation
	double	*Du; // divergence
	
	double	*p; // pressure
	double	*Gpx;
	double	*Gpy;
	double	*Gpz;
	double	*DGp;

	double	*bc_u_r; // velocity at right boundary
	double	*bc_v_r;
	double	*bc_w_r;
	double	*bc_u_l; // velocity at left boundary
	double	*bc_v_l;
	double	*bc_w_l;
	double	*bc_u_t; // velocity at top boundary
	double	*bc_v_t;
	double	*bc_w_t;
	double	*bc_u_b; // velocity at bottom boundary
	double	*bc_v_b;
	double	*bc_w_b;
	
	
	double *bc_dudy_t; 
	double *bc_rhs_dudy_t; // right hand side of wall ODE
	double *bc_dudy_b; 
	double *bc_rhs_dudy_b; // right hand side of wall ODE
	double *Hbc_dudy_t;
	double *Hbc_dudy_b;
	
	double *bc_kappa_t;
	double *bc_kappa_b;
	double *bc_K_t;
	double *bc_K_b;
	double *bc_txy_t;
	double *bc_txy_b;
	
	double *bc_convective_t;
	double *bc_convective_b;
	
	
	
	double *eta0_bar_t;
	double *eta0_t;
	double *lambda_t;
	double *eta0_bar_b;
	double *eta0_b;
	double *lambda_b;
	
	
	double *A1_mean;
	double *A2_mean;
	double *B1_mean;
	double *B2_mean;
	double *C1_mean;
	double *C2_mean;
	double *D_mean;
	double *E1_mean;
	double *E2_mean;
	

} variables;

/* ************************************************************************ */
/* Structure for naighboring data											*/
/* ************************************************************************ */

typedef struct neighbors{
	
	double *u_in_l;
	double *u_in_r;
	double *v_in_l;
	double *v_in_r;
	double *w_in_l;
	double *w_in_r;

} neighbors;

/* ************************************************************************ */
/* Structure for the fftw plans                                             */
/* ************************************************************************ */

typedef struct fftwplans{
	
	fftw_plan p1d_z_u;
	fftw_plan p1d_z_v;
	fftw_plan p1d_z_wp;
	fftw_plan p1d_z_bc_tb;
	fftw_plan p1d_z_bc_uw_lr;
	fftw_plan p1d_z_bc_v_lr;
	
	fftw_plan p1d_z_ued;
	fftw_plan p1d_z_ved;
	fftw_plan p1d_z_ct;
	fftw_plan p1d_z_cn;
	
	
	fftw_plan p1d_invz_u;
	fftw_plan p1d_invz_v;
	fftw_plan p1d_invz_wp;
	fftw_plan p1d_invz_bc_tb;
	fftw_plan p1d_invz_bc_uw_lr;
	fftw_plan p1d_invz_bc_v_lr;
	
	fftw_plan p1d_invz_ued;
	fftw_plan p1d_invz_ved;
	fftw_plan p1d_invz_ct;
	fftw_plan p1d_invz_cn;
	
	fftw_plan p1d_poisson;
	fftw_plan p1d_invpoisson;
	
	fftw_plan p2d_xz_wp; 
	fftw_plan p2d_invxz_wp; 
    
    fftw_plan transpose_yz2xy;
    fftw_plan transpose_xy2yz;
    
    fftw_plan p1d_x_wp;
    fftw_plan p1d_invx_wp;
	
    fftw_plan p1d_poisson_ver2;
    fftw_plan p1d_invpoisson_ver2;
} fftwplans;


/* ************************************************************************ */
/* Structure for interpolated values                                        */
/* ************************************************************************ */

typedef struct interpolated{
	
	double *u_ued;
	double *u_phys_ued;
	
	double *v_ved;
	double *v_phys_ved;
	
	double *w_ct;
	double *w_phys_ct;
	
	double *u_phys_les;
	double *v_phys_les;
	double *w_phys_les;
	
	
	double *u_phys_cn1;
	double *v_phys_cn1;
	
	double *u_phys_ct1;
	double *v_phys_ct1;
	
	double *w_phys_ved1;
	double *w_phys_ued1;
	
	double *u_phys_cn3;
	double *v_phys_cn3;
	
	double *u_phys_ct3;
	double *v_phys_ct3;
	
	double *w_phys_ved3;
	double *w_phys_ued3;
	
	
	double *uv3_top;
	double *uv3_bot;
	double *wv3_top;
	double *wv3_bot;
	
	double *bc_duudx_t;
	double *bc_duwdz_t;
	double *bc_duudx_b;
	double *bc_duwdz_b;
	
	double *bc_ududx_t;
	double *bc_wdudz_t;
	double *bc_ududx_b;
	double *bc_wdudz_b;
	
	double *bc_uv_t;
	double *bc_uv_b;
	
	
	
	double *Txx_ct;
	double *Tyy_ct;
	double *Tzz_ct;
	
	double *Txy_cn;
	double *Tyz_ved;
	double *Tzx_ued;
	
} interpolated;
