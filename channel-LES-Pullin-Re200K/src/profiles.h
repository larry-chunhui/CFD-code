extern void skin_friction(
                          double *u,
                          parameters *par );

extern void skin_friction__bc(
							  variables *var,
							  parameters *par );
extern void dump_bc_statistics(
                               variables *var,
                               parameters *par);

extern void statistics_u(
                         double *u,
                         double *v,
			 double *w,
                         parameters *par );

extern void mean_profile_eta(
							 double *res_diss,
							 parameters *par );

extern void mean_dissipation_ratio(
								   double *res_diss,
								   double *sub_diss,
								   parameters *par );
