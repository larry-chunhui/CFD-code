
	  IMPLICIT NONE

          INTEGER  E_ijk(3,2)
          
	  COMMON /E_ijk/E_ijk
	  
          REAL*8 ZL,XL,Re,DT

	  COMMON /MESH/ZL,XL,Re,DT

          REAL*8  G,QG1,Q



          REAL*8 time1,time2,TIME

          INTEGER  NT,it,IBC,ICONT,NSTEPS,IDUMFRQ,IPRNFRQ

	  INTEGER I,J,K,K1,K2,iy,iz,ix  

