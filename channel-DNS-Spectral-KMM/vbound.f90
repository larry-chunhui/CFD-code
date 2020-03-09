! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE VBOUND (VBCp,VBC,TIME)

	  USE FFT
	  INCLUDE 'dim.h'

	  REAL*8 VBCp(3,0:NX-1,0:1,0:NZ-1)
	  COMPLEX*16 VBC(3,0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1)

	  COMPLEX*16, ALLOCATABLE :: UK(:,:)
      ALLOCATE(UK(0:NX/2-1,-NZ/2+1:NZ/2-1))

	  VBCp = 0.0D0
!	  VBCp(1,NX/8,1,0) = 5.0D-2
!	  VBCp(2,NX/8,1,0) = 5.0D-2
!	  VBCp(1,7*NX/8,1,0) = 5.0D-2
!	  VBCp(2,7*NX/8,1,0) =-5.0D-2

	  DO I=1,3
		 DO J=0,1
			CALL xzfft (VBCp(I,0:NX-1,J,0:NZ-1),UK)
			VBC(I,0:NXP/2-1,J,-NZ/2+1:NZ/2-1)=UK(my_id*NXP/2:(my_id+1)*NXP/2-1,-NZ/2+1:NZ/2-1)
		 ENDDO
	  ENDDO

	  DEALLOCATE(UK)

      END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE xzfft (UP,U)

	  USE FFT

	  IMPLICIT REAL*8 (A-H, O-Z)

      REAL*8 UP(0:NX-1,0:NZ-1)
      COMPLEX*16 U(0:NX/2-1,-NZ/2+1:NZ/2-1)

	  REAL*8, ALLOCATABLE :: UMX(:),UK(:,:)
      COMPLEX*16, ALLOCATABLE :: UMZ(:)

      ALLOCATE(UMX(0:NX-1),UMZ(0:NZ-1),UK(0:NX-1,0:NZ-1))
	  
	  DO iz = 0,NZ-1
		 UMX = UP(0:NX-1,iz)
		 CALL rfftf(NX,UMX,trigx)
		 UK(0,iz) = UMX(0)/DFLOAT(NX)
		 UK(1,iz) = 0.0D0
		 UK(2:NX-1,iz) = UMX(1:NX-2)/DFLOAT(NX)
	  ENDDO

	  DO ix = 0,NX/2-1
		 DO iz = 0,NZ-1
			UMZ(iz) = DCMPLX(UK(2*ix,iz),UK(2*ix+1,iz))
		 ENDDO
		 CALL cfftf(NZ,UMZ,trigz)
		 DO iz = -NZ/2+1,-1
			U(ix,iz) = UMZ(NZ+iz)/DFLOAT(NZ)
		 ENDDO
		 DO iz = 0,NZ/2-1
			U(ix,iz) = UMZ(iz)/DFLOAT(NZ)
	     ENDDO
	  ENDDO

	  DEALLOCATE(UMX,UMZ)

      RETURN
	  END