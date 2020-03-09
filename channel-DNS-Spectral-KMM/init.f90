! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE INITIAL (Vp,VBCp)

	  USE FFT
	  INCLUDE 'dim.h'

	  REAL*8 Vp(3,0:NX-1,0:NY,0:NZP-1),VBCp(3,0:NX-1,0:1,0:NZ-1)
	  REAL*8 uy

	  Vp = 0.0D0

	  Vp(1:3,0:NX-1,0,0:NZP-1) = VBCp(1:3,0:NX-1,0,my_id*NZP:(my_id+1)*NZP-1)
      DO iy = 1,NY-1
!		 uy = yc(iy)
 		 uy = 1.5D0*(1.0D0-yc(iy)**2)
		 Vp(1,0:NX-1,iy,0:NZP-1) = uy
	  ENDDO
	  Vp(1:3,0:NX-1,NY,0:NZP-1) = VBCp(1:3,0:NX-1,1,my_id*NZP:(my_id+1)*NZP-1)

	  CALL DISTURB (Vp)

	  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE DISTURB (Vp)

	  USE FFT
	  INCLUDE 'dim.h'

	  REAL*8 Vp(3,0:NX-1,0:NY,0:NZP-1)
	  REAL*8 EPS,X,DRAND,uy
	  INTEGER NN

	  INTEGER DATE_TIME (8)
	  CHARACTER (LEN = 12) REAL_CLOCK (3)
	  CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), DATE_TIME)
	  !DATE_AND_TIME返回值有四组，DATE_TIME有8个值，其中7和8分别表示秒（0-60）和微秒（0-999）

	  EPS=0.1D0

	  CALL RANDOM_SEED

	  NN=100*DATE_TIME (8)+DATE_TIME (7)

	  X=DRAND(1)
	  DO iy=0,NN
		 X=DRAND(0)
	  ENDDO

      !在层流进口的基础上加一个随机扰动
      DO iy = 0,NY
		 uy = EPS*(1.0D0-yc(iy)**2)
		 DO iz = 0,NZP-1
			DO ix = 0,NX-1
			   Vp(1,ix,iy,iz) = Vp(1,ix,iy,iz)+uy*(DRAND(0)-0.5D0)
			   Vp(2,ix,iy,iz) = Vp(2,ix,iy,iz)+uy*(DRAND(0)-0.5D0)
			   Vp(3,ix,iy,iz) = Vp(3,ix,iy,iz)+uy*(DRAND(0)-0.5D0)
			ENDDO
		 ENDDO
	  ENDDO

	  END