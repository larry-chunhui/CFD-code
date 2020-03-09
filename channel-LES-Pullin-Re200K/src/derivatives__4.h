void get_divergence_bc__4(variables *var, parameters *par);

void get_laplace_bc_y__4(variables *var, parameters *par);
void get_laplace_bc_xz__4(variables *var, parameters *par);
void get_laplace_bc__4(variables *var, parameters *par);

void get_gradient__4(variables *var, int tag,
					 parameters *par);

void get_divergence__4(double *u, 
					   double *v, 
					   double *w, 
					   double *Du,
					   parameters *par);

void get_laplace_y__4(double *u,
					  double *Lu,
					  double *v,
					  double *Lv,
					  double *w,
					  double *Lw,
					  variables *var,
					  parameters *par);

void get_laplace_xz__4(double *u,
					   double *Lu,
					   double *v,
					   double *Lv,
					   double *w,
					   double *Lw,
					   variables *var,
					   parameters *par);

void get_laplace__4(double *u,
					double *Lu,
					double *v,
					double *Lv,
					double *w,
					double *Lw,
					variables *var,
					parameters *par);
