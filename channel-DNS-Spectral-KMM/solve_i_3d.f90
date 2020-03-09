! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE SOLVE_i (WN,AD,G,a,b)

	  USE FFT

	  IMPLICIT NONE

	  INTEGER ix,iy,iz
	  REAL*8 G,FC,a,b

      REAL*8 WN(0:NXP/2-1,-NZ/2+1:NZ/2-1),AD(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)

      DO iz = -NZ/2+1,NZ/2-1
         DO ix = 0,NXP/2-1
            WN(ix,iz) = -(WAVX(ix)**2 + WAVZ(iz)**2 + G)
         ENDDO
      ENDDO

      AD(0:NXP/2-1,ny,-NZ/2+1:NZ/2-1) = 1.0D0/(1.0D0+ WD(ny)*WN)
      AD(0:NXP/2-1,ny-1,-NZ/2+1:NZ/2-1) = 1.0D0/(1.0D0+ WD(ny-1)*WN)

      DO iy = NY-2,2,-1
         FC = WR(iy)*WL(iy+2)
         AD(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = 1.0D0/(1.0D0+(WD(iy)-FC*WN*AD(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1))*WN)
      ENDDO

	  RETURN
      END