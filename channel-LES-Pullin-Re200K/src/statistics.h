extern void statistics_les(
						   interpolated *inter,
						   fftwplans *planptr,
						   parameters *par );

extern void dump_reystress_les(
							   double *u_xz,
							   double *v_xz,
							   double *w_xz,
							   double *Txx_xz,
							   double *Tyy_xz,
							   double *Tzz_xz,
							   double *Txy_xz,
							   parameters *par);


extern void statistics(
                       interpolated *inter,
                       fftwplans *planptr,
                       parameters *par );

extern void dump_reystress(
                           double *u_xz,
                           double *v_xz,
                           double *w_xz,
                           parameters *par);


extern void profiles_reystress(
                               double *u_xz,
                               double *v_xz,
                               double *uv,
                               parameters *par );


void forward_xz(
                double *u_phys, 
                double *u_xz, 
                fftwplans *planptr,
                parameters *par);




extern void courant(
                    double *u_phys,
                    double *v_phys,
                    double *w_phys,
                    parameters *par );
