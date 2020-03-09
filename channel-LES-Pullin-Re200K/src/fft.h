#define FFTW_MALLOC( N ) ( ( double * )fftw_malloc( sizeof( double )*( N ) ) )

void fft_make_plans__transpose(parameters *par, fftwplans *planptr);
void fft_destroy_plans__transpose(fftwplans *planptr, parameters *par);
void fft_destroy_plans(fftwplans *planptr, parameters *par);
void fft_make_plans(variables *var, parameters *par, fftwplans *planptr);

extern void fftw_forward_poisson_ver2(double *phys, double *wave, fftwplans *planptr, parameters *par);
extern void fftw_inverse_poisson_ver2(double *phys, double *wave, fftwplans *planptr, parameters *par);
