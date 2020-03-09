extern void bc__time_march_bc_dudy_analitic(
											variables *var0,
											variables *var1,
											int substep,
											parameters *par );

extern void bc__time_march_bc_dudy(
								   variables *var_prev,
								   variables *var_curr,
								   variables *var_next,
								   parameters *par );
extern void bc__nonlinear(
						  variables *var,
						  interpolated *inter,
						  int substep,
						  parameters *par,
						  fftwplans *planptr);

extern void bc__rhs(
					variables *var_prev,
					variables *var_curr,
					int substep,
					parameters *par );


double average__nonlinear_ode(double *term, parameters *par);

extern void term_check__nonlinear_ode(
									  variables *var,
									  interpolated *inter,
									  int substep,
									  parameters *par,
									  fftwplans *planptr);

extern void wall_terms_check_dump(
								  double *A1, double *A2,
								  double *B1, double *B2,
								  double *C1, double *C2,
								  double *D,
								  double *E1, double *E2,
								  const char *basename,
								  parameters *par);

extern void mean_values_1d__term_check_nonlinear(
												 double *snap,
												 double *mean,
												 int substep,
												 parameters *par);


extern void bc_rhs_roughness(variables *var, double *ode_roughness_b, double *ode_roughness_t, parameters *par);
