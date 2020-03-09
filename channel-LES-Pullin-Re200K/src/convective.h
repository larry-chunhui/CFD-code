void trancate_exp(double *product, char *kinds, parameters *par);
void bc_trancate_exp(double *product, char *kinds, parameters *par);

void get_convective(
		variables *var,
		neighbors *shared,
		interpolated *inter,
		parameters *par,
		fftwplans *planptr);
		
void calc_multiply_ct(
		double *u_ct1,
		double *u_ct2, 
		double *uu,
		parameters *par);

void calc_multiply_ued(
		double *u_ed, 
		double *w_ued, 
		double *uw, 
		parameters *par);
		
void calc_multiply_ved(
		double *v_ed, 
		double *w_ved, 
		double *vw, 
		parameters *par);
		
void calc_multiply_cn(
		double *u_cn, 
		double *v_cn, 
		double *uv, 
		parameters *par);
