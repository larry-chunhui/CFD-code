! ---&---1---------2---------3---------4---------5---------6---------7--
!
	CHARACTER*2 FUNCTION itochar(iy)

	IMPLICIT REAL*8 (A-H, O-Z)

	WRITE(itochar,'(I2)') iy

	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	CHARACTER*8 FUNCTION Ftochar(Q)

	IMPLICIT REAL*8 (A-H, O-Z)

	WRITE(Ftochar,'(F8.2)') Q

	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE DOTTIME (ADB,A,B)

	USE FFT
	IMPLICIT REAL*8 (A-H,O-Z)

	COMPLEX*16 ADB(0:NX/2,0:NY-1,0:NZP-1)
	COMPLEX*16, DIMENSION (3,0:NX/2,0:NY-1,0:NZP-1) :: A,B

	REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ADB2(:,:,:)

	ALLOCATE (A2(3,0:NX2-1,0:NY2P-1,0:NZ2-1),B2(3,0:NX2-1,0:NY2P-1,0:NZ2-1),ADB2(0:NX2-1,0:NY2P-1,0:NZ2-1))

	DO I = 1,3
		CALL xyzfft2 ('B',A2(I,0:NX2-1,0:NY2P-1,0:NZ2-1),A(I,0:NX/2,0:NY-1,0:NZP-1))
		CALL xyzfft2 ('B',B2(I,0:NX2-1,0:NY2P-1,0:NZ2-1),B(I,0:NX/2,0:NY-1,0:NZP-1))
	ENDDO

	ADB2(0:NX2-1,0:NY2P-1,0:NZ2-1) =  A2(1,0:NX2-1,0:NY2P-1,0:NZ2-1)*B2(1,0:NX2-1,0:NY2P-1,0:NZ2-1) &
									+ A2(2,0:NX2-1,0:NY2P-1,0:NZ2-1)*B2(2,0:NX2-1,0:NY2P-1,0:NZ2-1) &
									+ A2(3,0:NX2-1,0:NY2P-1,0:NZ2-1)*B2(3,0:NX2-1,0:NY2P-1,0:NZ2-1)


	CALL xyzfft2 ('F',ADB2(0:NX2-1,0:NY2P-1,0:NZ2-1),ADB(0:NX/2,0:NY-1,0:NZP-1))

	DEALLOCATE (A2,B2,ADB2)

	RETURN
	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE CROSSTIME (ACB,A,B)

	USE FFT
	INCLUDE 'com_iso_3d.h'

	COMPLEX*16, DIMENSION (3,0:NX/2,0:NY-1,0:NZP-1) :: ACB,A,B

	REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ACB2(:,:,:,:)

	INTEGER I,K1,K2

	ALLOCATE (A2(3,0:NX2-1,0:NY2P-1,0:NZ2-1),B2(3,0:NX2-1,0:NY2P-1,0:NZ2-1),ACB2(3,0:NX2-1,0:NY2P-1,0:NZ2-1))



	DO I = 1,3
		CALL xyzfft2 ('B',A2(I,0:NX2-1,0:NY2P-1,0:NZ2-1),A(I,0:NX/2,0:NY-1,0:NZP-1))
		CALL xyzfft2 ('B',B2(I,0:NX2-1,0:NY2P-1,0:NZ2-1),B(I,0:NX/2,0:NY-1,0:NZP-1))
	ENDDO

	DO I = 1,3
		K1 = E_ijk(I,1);	K2 = E_ijk(I,2)
		ACB2(I,0:NX2-1,0:NY2P-1,0:NZ2-1) =	A2(K1,0:NX2-1,0:NY2P-1,0:NZ2-1)*B2(K2,0:NX2-1,0:NY2P-1,0:NZ2-1) - &
											A2(K2,0:NX2-1,0:NY2P-1,0:NZ2-1)*B2(K1,0:NX2-1,0:NY2P-1,0:NZ2-1)
	ENDDO

	DO I = 1,3
		CALL xyzfft2 ('F',ACB2(I,0:NX2-1,0:NY2P-1,0:NZ2-1),ACB(I,0:NX/2,0:NY-1,0:NZP-1))     
	ENDDO

	DEALLOCATE (A2,B2,ACB2)

	RETURN
	END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE Grad (Grad_U,U)

	USE FFT
	IMPLICIT REAL*8 (A-H,O-Z)

	COMPLEX*16 Grad_U(3,0:NX/2,0:NY-1,0:NZP-1),U(0:NX/2,0:NY-1,0:NZP-1)

	DO I = 1,3
		CALL DERIV (I,Grad_U(I,0:NX/2,0:NY-1,0:NZP-1),U(0:NX/2,0:NY-1,0:NZP-1))
	ENDDO

	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE Div(Div_U,U)

	USE FFT
	IMPLICIT REAL*8 (A-H,O-Z)

	COMPLEX*16 Div_U(0:NX/2,0:NY-1,0:NZP-1),U(3,0:NX/2,0:NY-1,0:NZP-1)

	COMPLEX*16, ALLOCATABLE :: DU(:,:,:,:)

	ALLOCATE (DU(3,0:NX/2,0:NY-1,0:NZP-1))

	DO I = 1,3
		CALL DERIV (I,DU(I,0:NX/2,0:NY-1,0:NZP-1),U(I,0:NX/2,0:NY-1,0:NZP-1))
	ENDDO

	Div_U = DU(1,0:NX/2,0:NY-1,0:NZP-1) + DU(2,0:NX/2,0:NY-1,0:NZP-1) + DU(3,0:NX/2,0:NY-1,0:NZP-1)

	DEALLOCATE (DU)

	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE Curl (Curl_U,U)

	USE FFT
	INCLUDE 'com_iso_3d.h'

	INTEGER I,K1,K2

	COMPLEX*16 Curl_U(3,0:NX/2,0:NY-1,0:NZP-1),U(3,0:NX/2,0:NY-1,0:NZP-1)

	COMPLEX*16, ALLOCATABLE :: DU(:,:,:,:)

	ALLOCATE (DU(2,0:NX/2,0:NY-1,0:NZP-1))

	DO I = 1,3
		K1 = E_ijk(I,1);	K2 = E_ijk(I,2)
		CALL DERIV (K1,DU(1,0:NX/2,0:NY-1,0:NZP-1),U(K2,0:NX/2,0:NY-1,0:NZP-1))
		CALL DERIV (K2,DU(2,0:NX/2,0:NY-1,0:NZP-1),U(K1,0:NX/2,0:NY-1,0:NZP-1))
		Curl_U(I,0:NX/2,0:NY-1,0:NZP-1) = DU(1,0:NX/2,0:NY-1,0:NZP-1) - DU(2,0:NX/2,0:NY-1,0:NZP-1)
	ENDDO

	DEALLOCATE (DU)

	END


! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE CROSSTIME2 (ACB,A,B)

	USE FFT
	INCLUDE 'com_iso_3d.h'

	COMPLEX*16, DIMENSION (3,0:NX/2,0:NY-1,0:NZP-1) :: ACB,A,B

	REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ACB2(:,:,:,:)

	INTEGER I,J,L,LL,K1,K2

	ALLOCATE (A2(3,0:NX-1,0:NYP-1,0:NZ-1),B2(3,0:NX-1,0:NYP-1,0:NZ-1),ACB2(3,0:NX-1,0:NYP-1,0:NZ-1))

!     write(*,*)"kmax=",KMAX
! 	WRITE(*,*)'WAVY=',WAVY(NY-KMAX-1)
! 	WRITE(*,*)'WAVZ=',WAVZ(NZ-KMAX-1)
! 	STOP
    DO I=KMAX+1,NX/2
	DO J=KMAX+1,NY-KMAX-1
	DO L=0,NZP-1
	   LL=my_id*NZP + L
	   IF(LL.GE.(KMAX+1).AND.LL.LE.(NZ-KMAX-1)) THEN
	      A(1:3,I,J,L)=DCMPLX(0.0D0,0.0D0)
	      B(1:3,I,J,L)=DCMPLX(0.0D0,0.0D0)
	   ENDIF
	ENDDO
	ENDDO
	ENDDO

! 	DO k = 0,NZP-1
! 		kk = my_id*NZP + k
! 		IF(kk>NZ/2) THEN
! 			WAVZ(k) = DFLOAT(kk-NZ)*ALPHAZ       !change 3.31
! 		ELSE
! 			WAVZ(k) = DFLOAT(kk)*ALPHAZ
! 		ENDIF
! 	ENDDO
    
	DO I = 1,3
		CALL xyzfft ('B',A2(I,0:NX-1,0:NYP-1,0:NZ-1),A(I,0:NX/2,0:NY-1,0:NZP-1))
		CALL xyzfft ('B',B2(I,0:NX-1,0:NYP-1,0:NZ-1),B(I,0:NX/2,0:NY-1,0:NZP-1))
	ENDDO

	DO I = 1,3
		K1 = E_ijk(I,1);	K2 = E_ijk(I,2)
		ACB2(I,0:NX-1,0:NYP-1,0:NZ-1) =	A2(K1,0:NX-1,0:NYP-1,0:NZ-1)*B2(K2,0:NX-1,0:NYP-1,0:NZ-1) - &
										A2(K2,0:NX-1,0:NYP-1,0:NZ-1)*B2(K1,0:NX-1,0:NYP-1,0:NZ-1)
	ENDDO

	DO I = 1,3
		CALL xyzfft ('F',ACB2(I,0:NX-1,0:NYP-1,0:NZ-1),ACB(I,0:NX/2,0:NY-1,0:NZP-1))     !”–Œ Ã‚
	ENDDO

	DEALLOCATE (A2,B2,ACB2)

	RETURN
	END

! ---&---1---------2---------3---------4---------5---------6---------7--