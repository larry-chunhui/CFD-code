void get_laplace_bc_y__2(variables *var, parameters *par);

void get_gradient__2(variables *var, int tag,
					 parameters *par);

void get_divergence__2(
					   double *u, 
					   double *v, 
					   double *w,
					   double *Du,
					   parameters *par);

void get_laplace_y__2(
					  double *u,
					  double *Lu,
					  double *v,
					  double *Lv,
					  double *w,
					  double *Lw,
					  variables *var,
					  parameters *par);

void get_laplace_xz__2(
					   double *u,
					   double *Lu,
					   double *v,
					   double *Lv,
					   double *w,
					   double *Lw,
					   variables *var,
					   parameters *par);
