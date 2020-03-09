! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE SETUP

	USE FFT
	INCLUDE 'com_iso_3d.h'
	REAL*8 ALPHAX,ALPHAY,ALPHAZ
	INTEGER k,kk,kx,ky,kz
	REAL*8 WAVEK,WAVEK2

	E_ijk(1,1) = 2;		E_ijk(1,2) = 3
	E_ijk(2,1) = 3;		E_ijk(2,2) = 1
	E_ijk(3,1) = 1;		E_ijk(3,2) = 2

	ALPHAX = 1.0D0/XL
	DO k = 0,NX/2
		WAVX(k) = DFLOAT(k)*ALPHAX
	ENDDO
	CWAVX = DCMPLX(0.0D0,1.0D0)*WAVX

	ALPHAY = 1.0D0/YL
	DO k = 0,NY-1                              !CHANGE 5.14
		IF(k>NY/2) THEN
			WAVY(k) = DFLOAT(k-NY)*ALPHAY       !change 3.31
		ELSE
			WAVY(k) = DFLOAT(k)*ALPHAY
		ENDIF
	ENDDO
	CWAVY = DCMPLX(0.0D0,1.0D0)*WAVY

	ALPHAZ = 1.0D0/ZL
	DO k = 0,NZP-1
		kk = my_id*NZP + k
		IF(kk>NZ/2) THEN
			WAVZ(k) = DFLOAT(kk-NZ)*ALPHAZ       !change 3.31
		ELSE
			WAVZ(k) = DFLOAT(kk)*ALPHAZ
		ENDIF
	ENDDO
	CWAVZ = DCMPLX(0.0D0,1.0D0)*WAVZ
	  
	DO kx = 0,NX/2
		DO ky = 0,NY-1
			DO kz = 0,NZP-1
				WAVEK2 = WAVX(kx)**2 + WAVY(ky)**2 + WAVZ(kz)**2
				IF(WAVEK2 < 1.0D-10) THEN
					EVK2(kx,ky,kz) = 1.0D0
					VECTORK(1:3,kx,ky,kz) = 0.0D0
				ELSE			
					EVK2(kx,ky,kz) = DEXP(DT*WAVEK2/Re)
					WAVEK = DSQRT(WAVEK2)
					VECTORK(1,kx,ky,kz) = WAVX(kx)/WAVEK
					VECTORK(2,kx,ky,kz) = WAVY(ky)/WAVEK
					VECTORK(3,kx,ky,kz) = WAVZ(kz)/WAVEK
				ENDIF
			ENDDO
		ENDDO
	ENDDO

	CALL rffti (NX,trigx)
	CALL cffti (NY,trigy)
	CALL cffti (NZ,trigz)
	CALL rffti (NX2,trigx2)
	CALL cffti (NY2,trigy2)
	CALL cffti (NZ2,trigz2)

	RETURN
	END