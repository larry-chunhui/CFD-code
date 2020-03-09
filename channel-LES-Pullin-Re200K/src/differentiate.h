void derivative_x1_ct2u(
			double *uu1,
			double *scratch, 
			parameters *par);

void derivative_x3_ct2u(
			double *uu3,
			double *scratch, 
			parameters *par);
	
void derivative_x1_cn2v(
			double *uv1,
			double *scratch, 
			parameters *par);
			
void derivative_x3_cn2v(
			double *uv3,
			double *scratch, 
			parameters *par);

void derivative_x1_ued2wp(
			double *uw1,
			double *scratch, 
			parameters *par);

void derivative_x3_ued2wp(
			double *uw3,
			double *scratch, 
			parameters *par);

void derivative_y1_cn2u(
			double *uv1,
			double *scratch, 
			parameters *par);
			
void derivative_y3_cn2u(
			double *uv3,
			double *scratch, 
			parameters *par);
			
void derivative_y1_ct2v(
			double *vv1,
			double *scratch, 
			parameters *par);
			
void derivative_y3_ct2v(
			double *vv3,
			double *scratch, 
			parameters *par);
	
void derivative_y1_ved2wp(
			double *vw1,
			double *scratch, 
			parameters *par);
			
void derivative_y3_ved2wp(
			double *vw3,
			double *scratch, 
			parameters *par);
	
void derivative_z_ued2u(
			double *uw,
			double *scratch, 
			parameters *par);

void derivative_z_ved2v(
			double *vw,
			double *scratch, 
			parameters *par);

void derivative_z_ct2wp(
			double *ww,
			double *scratch, 
			parameters *par);



void derivative_x1_ued2ct__2(
							 double *u,
							 double *dudx,
							 parameters *par);

void derivative_x1_ued2ct__4(
		double *u,
		double *dudx,
		parameters *par);
		
		
void derivative_x3_ued2ct__4(
		double *u,
		double *dudx,
		parameters *par);
		
void derivative_x1_ved2cn__2(
							 double *v,
							 double *dvdx,
							 parameters *par);

void derivative_x1_ved2cn__4(
		double *v,
		double *dvdx,
		parameters *par);
		
		
void derivative_x3_ved2cn__4(
		double *v,
		double *dvdx,
		parameters *par);
		
void derivative_x1_ct2ued__2(
							 double *w, 
							 double *dwdx, 
							 parameters *par);

void derivative_x1_ct2ued__4(
		double *w, 
		double *dwdx, 
		parameters *par);
		
		
void derivative_x3_ct2ued__4(
		double *w, 
		double *dwdx, 
		parameters *par);
		
void derivative_y1_ued2cn__2(
							 double *u, 
							 double *dudy,
							 parameters *par);

void derivative_y1_ued2cn__4(
		double *u, 
		double *dudy,
		parameters *par);
		
		
void derivative_y3_ued2cn__4(
		double *u, 
		double *dudy,
		parameters *par);

void derivative_y1_ved2ct__2(
							 double *v, 
							 double *dvdy, 
							 parameters *par);
		
void derivative_y1_ved2ct__4(
		double *v, 
		double *dvdy, 
		parameters *par);
		
		
void derivative_y3_ved2ct__4(
		double *v, 
		double *dvdy, 
		parameters *par);
	
void derivative_y1_ct2ved__2(
							 double *w, 
							 double *dwdy, 
							 parameters *par);
		
void derivative_y1_ct2ved__4(
		double *w, 
		double *dwdy, 
		parameters *par);
		
		
void derivative_y3_ct2ved__4(
		double *w, 
		double *dwdy, 
		parameters *par);
		
		
void derivative_z_ued2ued(
		double *u, 
		double *dudz, 
		parameters *par);
		
void derivative_z_ved2ved(
		double *v, 
		double *dvdz, 
		parameters *par);
		
void derivative_z_ct2ct(
		double *w,
		double *dwdz, 
		parameters *par);

void derivative_z_ued2ued__2(
						  double *u, 
						  double *dudz, 
						  parameters *par);

void derivative_z_ved2ved__2(
						  double *v, 
						  double *dvdz, 
						  parameters *par);

void derivative_z_ct2ct__2(
						double *w,
						double *dwdz, 
						parameters *par);
			
