extern void bc__get_bc_dudy(
							double *bc_u,
							double *bc_K,
							double *bc_dudy,
							char *kinds,
							parameters *par );


extern double bc__u_tau_from_u(
							   double u,
							   double gamma,
							   double K,
							   double z0p,
							   double nu,
							   double z );


extern void bc__get_bc_vel(
						   double *bc_dudz,
						   double *bc_K,
						   double *bc_kappa,
						   double *bc_u,
						   double *bc_v,
						   double *bc_w,
						   char *kinds,
						   parameters *par );


extern double bc__u_from_u_tau(
							   double u_tau,
							   double gamma,
							   double K,
							   double kappa,
							   double z0p,
							   double nu,
							   double z ,
                               parameters *par);
