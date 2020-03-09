! ---&---1---------2---------3---------4---------5---------6---------7--
!Ref.xinliang li page268
!CHANGE 15.4.7 CHUNHUI LIU
	SUBROUTINE INITIAL (V,A,K0)

	USE FFT
	INCLUDE 'com_iso_3d.h'

	COMPLEX*16 V(3,0:NX/2,0:NY-1,0:NZP-1)
	COMPLEX*16 U(3,0:NX/2,0:NY-1,0:NZP-1)


	INTEGER DATE_TIME (8)
	CHARACTER (LEN = 12) REAL_CLOCK (3)
	REAL*8 THETA      !(0,2PI)均匀分布的角度
	COMPLEX*16 COMLXTHETA,K_DOTU
	REAL*8 TWOPI,WAVEK2,DWAVEK2,XRAND,DRAND
	INTEGER kx,ky,kz
	REAL*8 EK

	REAL*8 K0,A

	CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), DATE_TIME)

	XRAND=DRAND(1000*DATE_TIME (8)+16*DATE_TIME (7))

	TWOPI = 2.0D0*DACOS(-1.0D0)
	  
!	AK02 = 4.0D0/DFLOAT(K0**2)

	DO kx = 0,NX/2
		DO ky = 0,NY-1
			DO kz = 0,NZP-1
				THETA= TWOPI*DRAND(0)
				COMLXTHETA = DCMPLX( DCOS(THETA) , DSIN(THETA) )
				IF(kx==0 .AND. ky==0.AND.kz==0) THEN
				V(1,kx,ky,kz)=DCMPLX(0.0D0,0.0D0)
				V(2,kx,ky,kz)=DCMPLX(0.0D0,0.0D0)
				V(3,kx,ky,kz)=DCMPLX(0.0D0,0.0D0)					
				ELSE
					WAVEK2 = WAVX(kx)**2 + WAVY(ky)**2 + WAVZ(kz)**2
					DWAVEK2 = DSQRT(WAVEK2)         
					CALL INI_EK(EK,A,K0,DWAVEK2)  
					U(1,kx,ky,kz)=DSQRT(EK/3.0D0/TWOPI/WAVEK2)*COMLXTHETA
					U(2,kx,ky,kz)=DSQRT(EK/3.0D0/TWOPI/WAVEK2)*COMLXTHETA
					U(3,kx,ky,kz)=DSQRT(EK/3.0D0/TWOPI/WAVEK2)*COMLXTHETA
					K_DOTU=VECTORK(1,kx,ky,kz)*U(1,kx,ky,kz)+&
					      +VECTORK(2,kx,ky,kz)*U(2,kx,ky,kz) +&
				          +VECTORK(3,kx,ky,kz)*U(3,kx,ky,kz)
				    V(1,kx,ky,kz)=VECTORK(1,kx,ky,kz)*K_DOTU-U(1,kx,ky,kz)
				    V(2,kx,ky,kz)=VECTORK(2,kx,ky,kz)*K_DOTU-U(2,kx,ky,kz)
				    V(3,kx,ky,kz)=VECTORK(3,kx,ky,kz)*K_DOTU-U(3,kx,ky,kz)
				ENDIF
			ENDDO
		ENDDO
	ENDDO

100 FORMAT(2F10.5,2F10.5,2F10.5)

	END

	SUBROUTINE INI_EK(E,A,K0,K)                              !ADD 3.24
	REAL*8 E,K,A,K0

!------ref li xinliang 
!	write(*,*) E,A,K0,K
	E=A*K**4*EXP(-2.0D0*(REAL(K)/REAL(K0))**2)
!	WRITE(*,*)E
!-------------------------

!------ref yuhaichuan
!    E=0.50D0*(REAL(K)**(-5.0D0/3.0D0))
	RETURN
	END