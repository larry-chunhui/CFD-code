void get_pentamatrix_uw(
						double *D, 
						int substep,
						parameters *par);

void get_pentamatrix_v(
					   double *D, 
					   int substep,
					   parameters *par);

void vel_pentasolver_u(
					   double *u_scratch, 
					   double *rhs_u, 
					   int substep,
					   parameters *par);

void vel_pentasolver_v(
					   double *v_scratch, 
					   double *rhs_v, 
					   int substep,
					   parameters *par);

void vel_pentasolver_w(
					   double *w_scratch, 
					   double *rhs_w, 
					   int substep,
					   parameters *par);
