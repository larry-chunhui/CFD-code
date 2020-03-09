! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE DERIV (IC,DF,F)

	  USE FFT

	  IMPLICIT NONE

      COMPLEX*16, DIMENSION(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) :: F,DF
	  INTEGER IC,ix,iy,iz

	  IF(IC==1) GO TO 10
	  IF(IC==2) GO TO 20
	  IF(IC==3) GO TO 30

  10  DO ix = 0,NXP/2-1
         DF(ix,0:NY,-NZ/2+1:NZ/2-1) = CWAVX(ix)*F(ix,0:NY,-NZ/2+1:NZ/2-1)
      ENDDO

      RETURN

 20   DF(0:NXP/2-1, ny ,-NZ/2+1:NZ/2-1) = DCMPLX(0.0D0, 0.0D0)
      DF(0:NXP/2-1,ny-1,-NZ/2+1:NZ/2-1) = DFLOAT(2*ny)*F(0:NXP/2-1,ny,-NZ/2+1:NZ/2-1)

      DO iy = NY-2,0,-1
         DF(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = DF(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1) &
									+ DFLOAT(2*iy+2)*F(0:NXP/2-1,iy+1,-NZ/2+1:NZ/2-1)
      ENDDO

      DF(0:NXP/2-1,0,-NZ/2+1:NZ/2-1) = DF(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)/2.0D0

      RETURN

 30   DO iz = -NZ/2+1,NZ/2-1
         DF(0:NXP/2-1,0:NY,iz) = CWAVZ(iz)*F(0:NXP/2-1,0:NY,iz)
      ENDDO

      RETURN
      END