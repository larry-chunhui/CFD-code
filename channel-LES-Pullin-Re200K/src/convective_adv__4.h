void get_convective_adv__4(
		interpolated *inter,
		double *Nu_copy,
		double *Nv_copy,
		double *Nw_copy,
		fftwplans *planptr,
		parameters *par);
		
void get_convective_adv__4_Nu(
		interpolated *inter,
		double *Nu_copy,
		fftwplans *planptr,
		parameters *par);
		
void get_convective_adv__4_Nv(
		interpolated *inter,
		double *Nv_copy,
		fftwplans *planptr,
		parameters *par);
		
void get_convective_adv__4_Nw(
		interpolated *inter,
		double *Nw_copy,
		fftwplans *planptr,
		parameters *par);
