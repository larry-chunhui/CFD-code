void get_convective_adv__2(
		interpolated *inter,
		double *Nu_copy,
		double *Nv_copy,
		double *Nw_copy,
		fftwplans *planptr,
		parameters *par);
		
void form_product_adi(
			interpolated *inter,
			double *ududx,
			double *udvdx,
			double *udwdx,
			double *vdudy,
			double *vdvdy,
			double *vdwdy,
			double *wdudz,
			double *wdvdz,
			double *wdwdz,
			fftwplans *planptr,
			parameters *par);
			
void form_product_x(
		interpolated *inter,
		double *ududx,
		double *udvdx,
		double *udwdx,
		fftwplans *planptr,
		parameters *par);

void form_product_y(
		interpolated *inter,
		double *vdudy,
		double *vdvdy,
		double *vdwdy,
		fftwplans *planptr,
		parameters *par);
		
void form_product_z(
		interpolated *inter,
		double *wdudz,
		double *wdvdz,
		double *wdwdz,
		fftwplans *planptr,
		parameters *par);

void get_dudx2(
		double *u,
		double *dudx,
		parameters *par);

void get_dvdx2(
		double *v,
		double *dvdx,
		parameters *par);

void get_dwdx2(
		double *w, 
		double *dwdx, 
		parameters *par);

void get_dudy2(
		double *u, 
		double *dudy,
		parameters *par);
		
void get_dvdy2(
		double *v, 
		double *dvdy, 
		parameters *par);

void get_dwdy2(
		double *w, 
		double *dwdy, 
		parameters *par);

void get_dudz_fft(
		double *u, 
		double *dudz, 
		parameters *par);
		
void get_dvdz_fft(
		double *v, 
		double *dvdz, 
		parameters *par);

void get_dwdz_fft(
		double *w, 
		double *dwdz, 
		parameters *par);


