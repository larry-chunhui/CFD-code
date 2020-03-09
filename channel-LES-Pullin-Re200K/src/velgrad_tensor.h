void get_velgrad_tensor_cn(
                           variables *var,
                           interpolated *inter,
                           parameters *par,
                           fftwplans *planptr);

void get_velgrad_tensor_cn__2(
							  variables *var,
							  interpolated *inter,
							  parameters *par,
							  fftwplans *planptr);


void get_dudy_cn_wall(double *u, double *u_t, double *u_b, double *dudy, parameters *par);
void get_dudy_wp_wall(double *u, double *u_t, double *u_b, double *dudy, parameters *par);
void get_dwdy_ved_wall(double *w, double *w_t, double *w_b, double *dwdy, parameters *par);
void get_dwdy_wp_wall(double *w, double *w_t, double *w_b, double *dwdy, parameters *par);
