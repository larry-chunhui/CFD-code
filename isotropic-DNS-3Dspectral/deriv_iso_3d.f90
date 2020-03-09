! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE DERIV (IC,DF,F)

	USE FFT
	IMPLICIT NONE
 
	COMPLEX*16, DIMENSION(0:NX/2,0:NY-1,0:NZP-1) :: F,DF
	INTEGER IC,ix,iy,iz

	IF(IC==1) GO TO 10
	IF(IC==2) GO TO 20
	IF(IC==3) GO TO 30

10	DO ix = 0,NX/2
		DF(ix,0:NY-1,0:NZP-1) = CWAVX(ix)*F(ix,0:NY-1,0:NZP-1)
	ENDDO

	RETURN

20	DO iy = 0,NY-1
		DF(0:NX/2,iy,0:NZP-1) = CWAVY(iy)*F(0:NX/2,iy,0:NZP-1)
	ENDDO

	RETURN

30	DO iz = 0,NZP-1
		DF(0:NX/2,0:NY-1,iz) = CWAVZ(iz)*F(0:NX/2,0:NY-1,iz)
	ENDDO

	RETURN
	END