void poisson_1DDCT_in_y_transpose_1DFFT_in_x(
                                             double *delta_p,
                                             double *rhs_b,
                                             parameters *par, 
                                             fftwplans *planptr);
			
void trisolver_helmholz_transpose_fft(
                                      double *delta_phat, 
                                      double *rhs_bhat, 
                                      parameters *par);


void septasolver_helmholz_transpose_fft(
                                        double *delta_phat, 
                                        double *rhs_bhat, 
                                        parameters *par);
					
