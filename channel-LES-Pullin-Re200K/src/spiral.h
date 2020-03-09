extern void spiral_sgs_stress(
                              double *u, double *v, double *w,
                              double *x, double *y, double *z,
                              int *yes_near, int size_near, int ix_this,
                              double dudx[3][3], double e[3], double nu, double del,
                              double *Txx, double *Tyy, double *Tzz,
                              double *Txy, double *Tyz, double *Tzx, double *prefac, double *lv);

extern void spiral_sgs_flux(
                            double dsdx[3], double e[3], double del, double K,
                            double *qx, double *qy, double *qz);

extern double spiral_ke_integral(double k);

extern double spiral_sf_integral(double d);

extern void spiral_eigenvalue_symm(
                                   double Sxx, double Syy, double Szz, double Sxy, double Syz, double Szx,
                                   double eigval[3]);

extern void spiral_eigenvector_symm(
                                    double Sxx, double Syy, double Szz, double Sxy, double Syz, double Szx,
                                    double eigval, double eigvec[3]);
