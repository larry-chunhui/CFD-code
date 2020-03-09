! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE SETUP

	  USE FFT
	  INCLUDE 'dim.h'   !°üº¬Re
	  REAL*8 PI,ALPHAX,ALPHAZ

	  E_ijk(1,1)=2;  E_ijk(1,2)=3
	  E_ijk(2,1)=3;  E_ijk(2,2)=1
	  E_ijk(3,1)=1;  E_ijk(3,2)=2

	  LSUB=1
	  LSUB(0)=0

	  L0=LSUB(my_id)

	  PI = 2.0D0*DACOS(0.0D0)

	  DO J = 0,NY/2-1
	     YC(J) = DCOS(PI*DFLOAT(J)/DFLOAT(NY))
	     YC(NY-J) = -YC(J)
	  ENDDO
	  YC(NY/2) = 0.0D0

      ALPHAX = 1.0D0/XL
      DO k = 0,NXP/2-1
         WAVX(k) = DFLOAT(my_id*NXP/2+k)*ALPHAX
      ENDDO
      CWAVX = DCMPLX(0.0D0,1.0D0)*WAVX
      CWAVXRe = CWAVX/Re

      ALPHAZ = 1.0D0/ZL
	  DO k = -NZ/2+1,NZ/2-1
         WAVZ(k) = DFLOAT(k)*ALPHAZ
      ENDDO
      CWAVZ = DCMPLX(0.0D0,1.0D0)*WAVZ
      CWAVZRe = CWAVZ/Re

	  CALL costi (NY+1,trigy)
	  CALL costi (NY2+1,trigy2)
	  CALL rffti (NX,trigx)
	  CALL rffti (NX2,trigx2)
	  CALL cffti (NZ,trigz)
	  CALL cffti (Nz2,trigz2)

 	  CALL rffti (NZ,trigzr)

      DO iy=NY,2,-1
         WL(iy) =  1.0D0/DFLOAT(4*iy*(iy-1))
         WD(iy) = -1.0D0/DFLOAT(2*(iy**2-1))
         WR(iy) =  1.0D0/DFLOAT(4*iy*(iy+1))
      ENDDO
      WL(2) = 2.0D0*WL(2)
      WD(ny-1:ny) = 0.0D0
      WR(ny-3:ny) = 0.0D0

	  RETURN
	  END