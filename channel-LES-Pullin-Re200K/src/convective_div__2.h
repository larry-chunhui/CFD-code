void get_convective_div__2(
		interpolated *inter,
		double *Nu_copy,
		double *Nv_copy,
		double *Nw_copy,
		fftwplans *planptr,
		parameters *par);
		

void form_product_div__2(
			interpolated *inter,
			double *uu,
			double *vv,
			double *ww,
			double *uv,
			double *uw,
			double *vw,
			parameters *par,
			fftwplans *planptr);
			
