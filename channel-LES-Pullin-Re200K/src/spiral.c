#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fasttrig.h"
#include "spiral.h"

#define EPS (2e-3)

#define EPS (2e-3)
#define TABLE_SIZE (int)1e5
#define DMAX 41.0
#define KMAX 40.0

static int initialized = 0;

static double spiral_ke_table[TABLE_SIZE];
static double spiral_sf_table[TABLE_SIZE];
void spiral_ke_LUT_init();
void spiral_sf_LUT_init();
double fastpow_spiral_sf_integral_ver3(double d);
double fastpow_spiral_ke_integral_ver3(double k);

//
// Calculate the SGS stresses (*Txx, *Tyy, *Tzz, *Txy, *Tyz, *Tzx) at
// (x[0], y[0], z[0]) given n (>=2) samples of the local resolved velocity
// field, (u[0], v[0], w[0]) at (x[0], y[0], z[0])
// to (u[n - 1], v[n - 1], w[n - 1]) at (x[n - 1], y[n - 1], z[n - 1]),
// resolved velocity gradient tensor dudx[3][3], LES cutoff scale del,
// and kinematic viscosity nu. If e[3] == { 0.0, 0.0, 0.0 }, overwrite e[3]
// with default alignment. prefac is the group prefactor
// \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3} and lv = Sqrt[2 nu / (3 Abs[a])],
// where a = e_i^v e_j^v S_{ij} is the axial stretching.
//
void spiral_sgs_stress(
                       double *u, double *v, double *w,
                       double *x, double *y, double *z,
                       int *yes_near, int size_near, int ix_this,
                       double dudx[3][3], double e[3], double nu, double del,
                       double *Txx, double *Tyy, double *Tzz,
                       double *Txy, double *Tyz, double *Tzx, double *K, double *lv)
{
    {
        // Strain-rate tensor
        double Sxx = 0.5 * (dudx[0][0] + dudx[0][0]);
        double Syy = 0.5 * (dudx[1][1] + dudx[1][1]);
        double Szz = 0.5 * (dudx[2][2] + dudx[2][2]);
        double Sxy = 0.5 * (dudx[0][1] + dudx[1][0]);
        double Syz = 0.5 * (dudx[1][2] + dudx[2][1]);
        double Szx = 0.5 * (dudx[2][0] + dudx[0][2]);
        
        if (e[0] == 0.0 && e[1] == 0.0 && e[2] == 0.0) {
            double eigval[3];
            spiral_eigenvalue_symm(
                                   Sxx, Syy, Szz, Sxy, Syz, Szx, eigval);
            // Default alignment: most extensive eigenvector
            spiral_eigenvector_symm(
                                    Sxx, Syy, Szz, Sxy, Syz, Szx, eigval[2], e);
		}
		
		
        // Make e[3] a unit vector
        double length = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        e[0] /= length; e[1] /= length; e[2] /= length;
        
        // Strain along vortex axis
        double a = e[0] * e[0] * Sxx + e[0] * e[1] * Sxy + e[0] * e[2] * Szx
        + e[1] * e[0] * Sxy + e[1] * e[1] * Syy + e[1] * e[2] * Syz
        + e[2] * e[0] * Szx + e[2] * e[1] * Syz + e[2] * e[2] * Szz;
        *lv = sqrt( 2.0 * nu / (3.0 * (fabs(a) + EPS)) );
    }
    
    if (initialized == 0){
        spiral_ke_LUT_init();
        spiral_sf_LUT_init();
        initialized = 1;
        printf("Initialized \n");
    }
    
    
    double F2 = 0.0;
    double Qd = 0.0;
    {
        // Average over neighboring points
        int i;
        for (i = 0; i < size_near; i++) {
			if ( yes_near[i] != 0) {
				double du = u[i] - u[ix_this];
				double dv = v[i] - v[ix_this];
				double dw = w[i] - w[ix_this];
				F2 += du * du + dv * dv + dw * dw;
				double dx = x[i] - x[ix_this];
				double dy = y[i] - y[ix_this];
				double dz = z[i] - z[ix_this];
				double dx2 = dx * dx   + dy * dy   + dz * dz;
				double dxe = dx * e[0] + dy * e[1] + dz * e[2];
				double d = sqrt(dx2 - dxe * dxe) / del;
				//Qd += spiral_sf_integral(d);
                Qd += fastpow_spiral_sf_integral_ver3(d); // faster version
			}
        }
        F2 /= (double) (size_near - 1);
        Qd /= (double) (size_near - 1);
    }
    double prefac = F2 / Qd; // \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
    double kc = M_PI / del;
    //*K = (prefac) * spiral_ke_integral(kc * (*lv));
    *K = (prefac) * fastpow_spiral_ke_integral_ver3(kc * (*lv)); // faster version
    
    // T_{ij} = (\delta_{ij} - e_i^v e_j^v) K
    *Txx = (1.0 - e[0] * e[0]) * (*K);
    *Tyy = (1.0 - e[1] * e[1]) * (*K);
    *Tzz = (1.0 - e[2] * e[2]) * (*K);
    *Txy = (    - e[0] * e[1]) * (*K);
    *Tyz = (    - e[1] * e[2]) * (*K);
    *Tzx = (    - e[2] * e[0]) * (*K);
}

//
// Calculate the SGS scalar flux (*qx, *qy, *qz) given the resolved scalar
// gradient dsdx[3], vortex alignment e[3] (a unit vector), LES cutoff
// scale del and precalculated SGS kinetic energy K.
//
//void spiral_sgs_flux(
//    double dsdx[3], double e[3], double del, double K,
//    double *qx, double *qy, double *qz)
//{
//    double gam = 1.0; // Universal model constant
//    double P = -0.5 * gam * del * sqrt(K);
//
//    // q_i = P (\delta_{ij} - e_i^v e_j^v) ds/dx_j
//    *qx = P * ((1.0 - e[0] * e[0]) * dsdx[0]
//             + (    - e[0] * e[1]) * dsdx[1]
//             + (    - e[0] * e[2]) * dsdx[2]);
//    *qy = P * ((    - e[1] * e[0]) * dsdx[0]
//             + (1.0 - e[1] * e[1]) * dsdx[1]
//             + (    - e[1] * e[2]) * dsdx[2]);
//    *qz = P * ((    - e[2] * e[0]) * dsdx[0]
//             + (    - e[2] * e[1]) * dsdx[1]
//             + (1.0 - e[2] * e[2]) * dsdx[2]);
//}

//
// Approximation of
// (1/2) k^(2/3) Gamma[-1/3, k^2]
// with maximum relative error of 0.17% at k=2.42806.
//
double spiral_ke_integral(double k)
{
    double k2 = k * k;
    if (k2 < 2.42806) {
        double pade = (3.0 +   2.5107 * k2 +  0.330357 * k2 * k2
                       +  0.0295481 * k2 * k2 * k2)
        / (1.0 + 0.336901 * k2 + 0.0416684 * k2 * k2
           + 0.00187191 * k2 * k2 * k2);
        return 0.5 * (pade - 4.06235 * pow(k2, 1.0 / 3.0));
    }
    else {
        double pade = (1.26429 + 0.835714 * k2 + 0.0964286 * k2 * k2)
        / (1.0     +   2.25   * k2 +  0.964286 * k2 * k2
           + 0.0964286 * k2 * k2 * k2);
        return 0.5 * pade * exp(-k2);
    }
}

//
// Approximation of
// Integrate[4 x^(-5/3) (1 - BesselJ[0, x Pi d]), {x, 0, 1}]
// with maximum relative error of 2.71% at d=0.873469.
//
double spiral_sf_integral(double d)
{
    // Uncomment if spherical averaging and d=1.
    // if (d == 1.0) return 4.09047;
    
    double d2 = d * d;
    if (d < 0.873469)
        return 7.4022 * d2 - 1.82642 * d2 * d2;
    else
        return 12.2946 * pow(d, 2.0 / 3.0) - 6.0
        - 0.573159 * pow(d, -1.5) * sin(3.14159 * d - 0.785398);
}

//
// Calculate the eigenvalues, eigval[0] < eigval[1] < eigval[2],
// of the 3 x 3 symmetric matrix,
// { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
// assuming distinct eigenvalues.
//
void spiral_eigenvalue_symm(
                            double Sxx, double Syy, double Szz, double Sxy, double Syz, double Szx,
                            double eigval[3])
{
    // x^3 + a * x^2 + b * x + c = 0, where x is the eigenvalue
    double a = - (Sxx + Syy + Szz);
    double b = Sxx * Syy - Sxy * Sxy + Syy * Szz
    - Syz * Syz + Szz * Sxx - Szx * Szx;
    double c = - (Sxx * (Syy * Szz - Syz * Syz)
                  + Sxy * (Syz * Szx - Sxy * Szz)
                  + Szx * (Sxy * Syz - Syy * Szx));
    
    double q = (3.0 * b - a * a) / 9.0;
    double r = (9.0 * a * b - 27.0 * c - 2.0 * a * a * a) / 54.0;
    
    if (q >= 0.0) {
        fprintf(stderr, "spiral_eigenvalue_symm(): q >= 0.0\n");
        exit(EXIT_FAILURE);
    }
    
    double costheta = r / sqrt(-q * q * q);
    
    // |costheta| > 1 should not occur, except from round-off errors
    double theta = costheta > 1.0 ? theta = 0.0 :
    costheta < -1.0 ? theta = M_PI :
    acos(costheta);
    
    eigval[0] = 2.0 * sqrt(-q) * cos((theta             ) / 3.0) - a / 3.0;
    eigval[1] = 2.0 * sqrt(-q) * cos((theta + 2.0 * M_PI) / 3.0) - a / 3.0;
    eigval[2] = 2.0 * sqrt(-q) * cos((theta + 4.0 * M_PI) / 3.0) - a / 3.0;
    
    // Sort eigenvalues: eigval[0] < eigval[1] < eigval[2]
    if (eigval[0] > eigval[1]) {
        double tmp = eigval[0]; eigval[0] = eigval[1]; eigval[1] = tmp;
    }
    if (eigval[1] > eigval[2]) {
        double tmp = eigval[1]; eigval[1] = eigval[2]; eigval[2] = tmp;
    }
    if (eigval[0] > eigval[1]) {
        double tmp = eigval[0]; eigval[0] = eigval[1]; eigval[1] = tmp;
    }
}

//
// Calculate the eigenvector (not normalized), eigvec[3],
// corresponding to the precalculated eigenvalue, eigval,
// of the 3 x 3 symmetric matrix,
// { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
// assuming distinct eigenvalues.
//
void spiral_eigenvector_symm(
                             double Sxx, double Syy, double Szz, double Sxy, double Syz, double Szx,
                             double eigval, double eigvec[3])
{
	
	// normalized by Frobenius norm
	double norm = sqrt(Sxx*Sxx+Syy*Syy+Szz*Szz+Sxy*Sxy+Syz*Syz+Szx*Szx);
	
	if (fabs((Sxx - eigval) * ((Syy - eigval) * (Szz - eigval) - Syz * Syz)
             + Sxy * (Syz * Szx - Sxy * (Szz - eigval))
             + Szx * (Sxy * Syz - (Syy - eigval) * Szx))/fabs(norm) > EPS) {
        fprintf(stderr, "spiral_eigenvector_symm(): invalid eigenvalue\n");
        exit(EXIT_FAILURE);
    }
	
    double det[3] = { (Syy - eigval) * (Szz - eigval) - Syz * Syz,
        (Szz - eigval) * (Sxx - eigval) - Szx * Szx,
        (Sxx - eigval) * (Syy - eigval) - Sxy * Sxy };
    
    double fabsdet[3] = { fabs(det[0]), fabs(det[1]), fabs(det[2]) };
    
    if (fabsdet[0] >= fabsdet[1] && fabsdet[0] >= fabsdet[2]) {
        eigvec[0] = 1.0;
        eigvec[1] = (-Sxy * (Szz - eigval) + Szx * Syz) / det[0];
        eigvec[2] = (-Szx * (Syy - eigval) + Sxy * Syz) / det[0];
    }
    else if (fabsdet[1] >= fabsdet[2] && fabsdet[1] >= fabsdet[0]) {
        eigvec[0] = (-Sxy * (Szz - eigval) + Syz * Szx) / det[1];
        eigvec[1] = 1.0;
        eigvec[2] = (-Syz * (Sxx - eigval) + Sxy * Szx) / det[1];
    }
    else if (fabsdet[2] >= fabsdet[0] && fabsdet[2] >= fabsdet[1]) {
        eigvec[0] = (-Szx * (Syy - eigval) + Syz * Sxy) / det[2];
        eigvec[1] = (-Syz * (Sxx - eigval) + Szx * Sxy) / det[2];
        eigvec[2] = 1.0;
    }
    else {
        fprintf(stderr, "spiral_eigenvector_symm():\n");
        exit(EXIT_FAILURE);
    }
}





void spiral_ke_LUT_init()
{
    long int i;
    for (i=0; i<TABLE_SIZE; i++) {
        double k = (KMAX*i/TABLE_SIZE);
        double k2 = k * k;
        double k4 = k2 * k2;
        
        if (k2 < 2.42806) {
            double pade = (3.0 +   2.5107 * k2 +  0.330357 * k4
                           +  0.0295481 * k2 * k4)
            / (1.0 + 0.336901 * k2 + 0.0416684 * k4
               + 0.00187191 * k2 * k4);
            spiral_ke_table[i] = 0.5 * (pade - 4.06235 * pow(k2, 1.0 / 3.0));
        }
        else {
            double pade = (1.26429 + 0.835714 * k2 + 0.0964286 * k4)
            / (1.0     +   2.25   * k2 +  0.964286 * k4
               + 0.0964286 * k2 * k4);
            spiral_ke_table[i] = 0.5 * pade * exp(-k2);
        }
    }
}



void spiral_sf_LUT_init()
{
    long int i;
    for (i=0; i<TABLE_SIZE; i++) {
        double d = (DMAX*i/TABLE_SIZE);
        
        
        double theta = 3.14159 * d - 0.785398;
        
        double d2 = d * d;
        
        if (d < 0.873469)
            spiral_sf_table[i] =  7.4022 * d2 - 1.82642 * d2 * d2;
        else
            spiral_sf_table[i] = 12.2946 * pow(d, 2.0 / 3.0) - 6.0
            - 0.573159 * pow(d, -1.5) * sin(theta);
        
    }
}



double fastpow_spiral_sf_integral_ver3(double d)
{
    
    if (d > DMAX){
        printf("fastpow_spiral_sf_integral_ver3 d = %e\n", d);
        return 12.2946 * pow(d, 2.0 / 3.0) - 6.0
        - 0.573159 * pow(d, -1.5) * sin(3.14159 * d - 0.785398); 
        
    }else if (d < 0.873469){
        double d2 = d * d;
        return 7.4022 * d2 - 1.82642 * d2 * d2;
    }else{
        long int pow_idx = (TABLE_SIZE/DMAX)*d;
        
        double d0 = DMAX/TABLE_SIZE*pow_idx;
        double d1 = DMAX/TABLE_SIZE*(pow_idx+1);
        
        double sfd0 = spiral_sf_table[pow_idx];
        double sfd1 = spiral_sf_table[pow_idx+1];
        
        
        return sfd0 + (d-d0)*(sfd1-sfd0)/(d1-d0);
    }
}



double fastpow_spiral_ke_integral_ver3(double k)
{
    
    if (k > KMAX){
        
        
        double k2 = k * k;
        double k4 = k2 * k2;
        
        printf("fastpow_spiral_ke_integral_ver3 k2 = %e\n", k2);
        double pade = (1.26429 + 0.835714 * k2 + 0.0964286 * k4)
        / (1.0     +   2.25   * k2 +  0.964286 * k4
           + 0.0964286 * k2 * k4);
        return 0.5 * pade * exp(-k2);
        
    }else{
        
        long int pow_idx = (TABLE_SIZE/KMAX)*k;
        
        double k0 = KMAX/TABLE_SIZE*pow_idx;
        double k1 = KMAX/TABLE_SIZE*(pow_idx+1);
        double kek0 = spiral_ke_table[pow_idx];
        double kek1 = spiral_ke_table[pow_idx+1];
        
        return kek0 + (k-k0)*(kek1-kek0)/(k1-k0);
        
    }
}
