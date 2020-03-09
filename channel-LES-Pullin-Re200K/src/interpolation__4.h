void interpolation__4(
		int flag,
		variables *var,
		neighbors *shared,
		interpolated *inter, 
		parameters *par,
		fftwplans *planptr);

		
void interpolate_y1_ct2ved__4(
		double *w_ct,
		double *w_ved,
		parameters *par);
		
void interpolate_x1_ct2ued__4(
		double *w_ct,
		double *w_ued,
		parameters *par);
		
void interpolate_x1_ved2cn__4(
		double *v_ved,
		double *v_cn,
		parameters *par);
		
void interpolate_y1_ved2ct__4(
		double *v_ved,
		double *v_ct,
		parameters *par);
		
void interpolate_y1_ued2cn__4(
		double *u_ued,
		double *u_cn,
		parameters *par);

void interpolate_x1_ued2ct__4(
		double *u_ued,
		double *u_ct,
		parameters *par);
		
void interpolate_y3_ct2ved__4(
		double *w_ct,
		double *w_ved,
		parameters *par);
		
void interpolate_x3_ct2ued__4(
		double *w_ct,
		double *w_ued,
		parameters *par);
		
void interpolate_x3_ved2cn__4(
		double *v_ved,
		double *v_cn,
		parameters *par);
		
void interpolate_y3_ved2ct__4(
		double *v_ved,
		double *v_ct,
		parameters *par);
		
void interpolate_y3_ued2cn__4(
		double *u_ued,
		double *u_cn,
		parameters *par);

void interpolate_x3_ued2ct__4(
		double *u_ued,
		double *u_ct,
		parameters *par);
		
		


void interpolate_y1_ct2v__4(
		double *w_ct,
		double *w_v,
		parameters *par);
		
void interpolate_x1_ct2u__4(
		double *w_ct,
		double *w_u,
		parameters *par);
		
void interpolate_y1_ved2wp__4(
		double *v_ved,
		double *v_wp,
		parameters *par);

void interpolate_x1_ued2wp__4(
		double *u_ued,
		double *u_wp,
		parameters *par);
		
void interpolate_y3_ct2v__4(
		double *w_ct,
		double *w_v,
		parameters *par);
		
void interpolate_x3_ct2u__4(
		double *w_ct,
		double *w_u,
		parameters *par);
		
void interpolate_y3_ved2wp__4(
		double *v_ved,
		double *v_wp,
		parameters *par);

void interpolate_x3_ued2wp__4(
		double *u_ued,
		double *u_wp,
		parameters *par);



void interpolate_y1_cn2u__4(
		double *v_cn,
		double *v,
		parameters *par);
	
void interpolate_y3_cn2u__4(
		double *v_cn,
		double *v,
		parameters *par);


void interpolate_y1_cn2ued__4(
		double *u_cn,
		double *u_ed,
		parameters *par);

void interpolate_y3_cn2ued__4(
		double *u_cn,
		double *u_ed,
		parameters *par);


void interpolate_x1_cn2v__4(
		double *v_cn,
		double *v,
		parameters *par);

void interpolate_x3_cn2v__4(
		double *v_cn,
		double *v,
		parameters *par);

void interpolate_x1_cn2ved__4(
		double *v_cn,
		double *ved,
		parameters *par);

void interpolate_x3_cn2ved__4(
		double *v_cn,
		double *ved,
		parameters *par);



void bigger_array_u_ued__4(
		double flag,
		variables *var,
		neighbors *shared,
		double *u_ed,
		parameters *par);

void bigger_array_v_ved__4(
		double flag,
		variables *var,
		neighbors *shared,
		double *v_ed,
		parameters *par);
		
void bigger_array_wp_ct__4(
		double flag,
		variables *var,
		neighbors *shared,
		double *w_ct,
		parameters *par);

void bigger_array_v_ved__4_addition(
					double *v_ed,
					double *u_ed,
					double *w_ct,
					parameters *par);

void interpolate_ct2cn__4(
						  double *ct,
						  double *cn,
						  parameters *par);

void interpolate_cn2wp__4(
						  double *cn,
						  double *wp,
						  parameters *par);

