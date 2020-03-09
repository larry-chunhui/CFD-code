void get_rhs_u(
			variables *var1,
			variables *var0,
			variables *var9,
			double *rhs_u,
			int substep,
			parameters *par);
					
void get_rhs_v(
			variables *var1,
			variables *var0,
			variables *var9,
			double *rhs_v,
			int substep,
			parameters *par);
					
void get_rhs_w(
			variables *var1,
			variables *var0,
			variables *var9,
			double *rhs_w,
			int substep,
			parameters *par);


void get_rhs_poisson(
			variables *var1,
			double *rhs_poisson,
			int substep,
			parameters *par);
