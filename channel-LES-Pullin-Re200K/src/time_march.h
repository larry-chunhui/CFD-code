void time_march(
				variables *var9,
				variables *var0,
				variables *var1,
				parameters *par,
				int substep,
				fftwplans *planptr);
				
void get_rhs(
		variables *var1,
		variables *var0,
		variables *var9,
		double *rhs_u,
		double *rhs_v,
		double *rhs_w,
		int substep,
		parameters *par);
		
void vel_solvers(
		double *u,
		double *v,
		double *w,
		double *rhs_u,
		double *rhs_v,
		double *rhs_w,
		int substep,
		parameters *par);
		

void get_rhs_transpose(
		variables *var1,
		variables *var0,
		variables *var9,
		double *rhs_u,
		double *rhs_v,
		double *rhs_w,
		int substep,
		parameters *par);
		
void vel_solvers_transpose(
		double *u,
		double *v,
		double *w,
		double *rhs_u,
		double *rhs_v,
		double *rhs_w,
		int substep,
		parameters *par);


void get_dTijdxj__4(
				 variables *var,
				 interpolated *inter,
				 parameters *par,
				 fftwplans *planptr);

void get_dTijdxj__2(
					variables *var,
					interpolated *inter,
					parameters *par,
					fftwplans *planptr);

void get_dTijdxj__test(
					   variables *var,
					   interpolated *inter,
					   double *dTxxdx,
					   double *dTxydy,
					   double *dTzxdz,
					   double *dTxydx,
					   double *dTyydy,
					   double *dTyzdz,
					   double *dTzxdx,
					   double *dTyzdy,
					   double *dTzzdz,
					   parameters *par,
					   fftwplans *planptr);
void get_dTijdxj__2_test(
					   variables *var,
					   interpolated *inter,
					   double *dTxxdx,
					   double *dTxydy,
					   double *dTzxdz,
					   double *dTxydx,
					   double *dTyydy,
					   double *dTyzdz,
					   double *dTzxdx,
					   double *dTyzdy,
					   double *dTzzdz,
					   parameters *par,
					   fftwplans *planptr);

void set_Tij_zero(variables *var, interpolated *inter, parameters *par);
