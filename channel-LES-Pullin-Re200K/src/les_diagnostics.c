/*
 *  les_diagnostics.c
 *  3D channel-flow, periodic in x-direction version.
 *	non-parallel version
 *
 */

#include "definitions.h"
#include "les_diagnostics.h"

extern void les_get_dissipation(
                                double *dudx, double *dudy, double *dudz,
                                double *dvdx, double *dvdy, double *dvdz,
                                double *dwdx, double *dwdy, double *dwdz,
                                double *Txx,  double *Tyy,  double *Tzz,
                                double *Txy,  double *Tyz,  double *Tzx,
                                double *res_diss,
                                double *sub_diss,
                                parameters *par )
{
    
	int i, j, k;
	int index;
	
	const int ny = par->ny;
	const int nz = par->nz;
	const int start = par->local_nx_start;
	const int end = par->local_nx_start + par->local_nx;
	
    double duidxj[3][3], sij[3][3], tij[3][3], nu, res_eps=0, sub_eps=0;
    int p, q;
    
    
    nu = 1.0/par->re;
    
    for (i=start; i<=end-1; ++i) {
        for (j=0; j<=ny; ++j) {
            for (k=0; k<=nz-1; ++k) {
                
				index = par->index_cn[i][j][k];
                
                duidxj[0][0] = dudx[index];
                duidxj[0][1] = dudy[index];
                duidxj[0][2] = dudz[index];
                duidxj[1][0] = dvdx[index];
                duidxj[1][1] = dvdy[index];
                duidxj[1][2] = dvdz[index];
                duidxj[2][0] = dwdx[index];
                duidxj[2][1] = dwdy[index];
                duidxj[2][2] = dwdz[index];
                
                tij[0][0] = Txx[index];
                tij[0][1] = Txy[index];
                tij[0][2] = Tzx[index];
                tij[1][0] = Txy[index];
                tij[1][1] = Tyy[index];
                tij[1][2] = Tyz[index];
                tij[2][0] = Tzx[index];
                tij[2][1] = Tyz[index];
                tij[2][2] = Tzz[index];
                
                for ( p = 0; p < 3; p++ ) {
                    for ( q = 0; q < 3; q++ ) {
                        sij[p][q] = ( duidxj[p][q] + duidxj[q][p] )/2.0;
					}
				}
                
                res_eps = 0.0;
                for ( p = 0; p < 3; p++ ) {
                    for ( q = 0; q < 3; q++ ) {
                        res_eps += sij[p][q]*sij[p][q];
					}
				}
                res_eps *= 2.0*nu;
                
                for ( p = 0; p < 3; p++ ) {
                    for ( q = 0; q < 3; q++ ) {
                        sub_eps += sij[p][q]*tij[p][q];
					}
				}
                
                sub_eps *= -1.0;
				
                
                res_diss[index] = res_eps;
                sub_diss[index] = sub_eps;
            }
		}
	}
	
	return;
	
}
