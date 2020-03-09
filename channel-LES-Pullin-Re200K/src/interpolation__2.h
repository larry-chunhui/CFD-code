void interpolation__2(
					  int flag,
					  variables *var,
					  neighbors *shared,
					  interpolated *inter, 
					  parameters *par,
					  fftwplans *planptr);

void interpolate1_ct_ved__2(
							double *w_ct,
							double *w_ved,
							parameters *par);

void interpolate1_ct_ued__2(
							double *w_ct,
							double *w_ued,
							parameters *par);

void interpolate1_ved_cn__2(
							double *v_ved,
							double *v_cn,
							parameters *par);

void interpolate1_ved_ct__2(
							double *v_ved,
							double *v_ct,
							parameters *par);

void interpolate1_ued_cn__2(
							double *u_ued,
							double *u_cn,
							parameters *par);

void interpolate1_ued_ct__2(
							double *u_ued,
							double *u_ct,
							parameters *par);

void interpolate_y1_cn2ued__2(
							  double *u_cn,
							  double *u_ed,
							  parameters *par);

void interpolate_x1_cn2ved__2(
							  double *v_cn,
							  double *ved,
							  parameters *par);

void interpolate_x1_ued2wp__2(
							  double *u_ued,
							  double *u_wp,
							  parameters *par);

void interpolate_cn2wp__2(
						  double *cn,
						  double *wp,
						  parameters *par);

void interpolate_ct2cn__2(
						  double *ct,
						  double *cn,
						  parameters *par);

void bigger_array_u_ued__2(
						   double flag,
						   variables *var,
						   neighbors *shared,
						   double *u_ed,
						   parameters *par);

void bigger_array_v_ved__2(
						   double flag,
						   variables *var,
						   neighbors *shared,
						   double *v_ed,
						   parameters *par);

void bigger_array_wp_ct__2(
						   double flag,
						   variables *var,
						   neighbors *shared,
						   double *w_ct,
						   parameters *par);

