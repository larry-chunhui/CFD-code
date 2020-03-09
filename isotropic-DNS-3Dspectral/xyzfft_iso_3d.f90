! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE xyzfft (C,UP,U)

	USE FFT
	IMPLICIT NONE

	include 'mpif.h'
	integer status(MPI_Status_Size)
	integer ierr

	CHARACTER C
	REAL*8 UP(0:NX-1,0:NYP-1,0:NZ-1)
	COMPLEX*16 U(0:NX/2,0:NY-1,0:NZP-1)
	INTEGER iy,iz,ix,K

	REAL*8, ALLOCATABLE :: UMX(:)
	COMPLEX*16, ALLOCATABLE :: UMY(:),UMZ(:),UK(:,:,:),UY(:,:,:)

	ALLOCATE (UMX(0:NX-1),UMY(0:NY-1),UMZ(0:NZ-1),UK(0:NX/2,0:NYP-1,0:NZ-1),UY(0:NX/2,0:NY-1,0:NZP-1))


	IF(C.EQ.'F') GO TO 10
	IF(C.EQ.'B') GO TO 20

! ---
10	DO iy = 0,NYP-1
		DO iz = 0,NZ-1
			UMX = UP(0:NX-1,iy,iz)
			CALL rfftf(NX,UMX,trigx)
			UK(0,iy,iz) = DCMPLX( UMX(0)/DFLOAT(NX) , 0.0D0 )
			DO ix = 1,NX/2-1
				UK(ix,iy,iz) = DCMPLX( UMX(2*ix-1)/DFLOAT(NX) , UMX(2*ix)/DFLOAT(NX) )
			ENDDO
			IF(INDEX_NXYZ==0) THEN
				UK(NX/2,iy,iz) = DCMPLX( 0.0D0 , 0.0D0 )
			ELSE
				UK(NX/2,iy,iz) = DCMPLX( UMX(NX-1)/DFLOAT(NX+NX) , 0.0D0 )
			ENDIF
		ENDDO
	ENDDO
	  
	DO ix = 0,NX/2
		DO iy = 0,NYP-1
			UMZ = UK(ix,iy,0:NZ-1)
			CALL cfftf(NZ,UMZ,trigz)
			UK(ix,iy,0:NZ/2-1) = UMZ(0:NZ/2-1)/DFLOAT(NZ)
			IF(INDEX_NXYZ==0) THEN
				UK(ix,iy,NZ/2) = DCMPLX( 0.0D0 , 0.0D0 )
			ELSE
				UK(ix,iy,NZ/2) = UMZ(NZ/2)/DFLOAT(NZ+NZ)
			ENDIF
			UK(ix,iy,NZ/2+1:NZ-1) = UMZ(NZ/2+1:NZ-1)/DFLOAT(NZ)
		ENDDO
	ENDDO

	IF(INDEX_Serial==1) THEN
		UY = UK
	ELSE
		DO K = 0,NPROC-1
			CALL MPI_SENDRECV(UK(0:NX/2,0:NYP-1,K*NZP:(K+1)*NZP-1),(NX/2+1)*NYP*NZP,MPI_COMPLEX16,K,10, &
							UY(0:NX/2,K*NYP:(K+1)*NYP-1,0:NZP-1),(NX/2+1)*NYP*NZP,MPI_COMPLEX16,K,10,MPI_COMM_WORLD,STATUS,ierr)
		ENDDO

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	ENDIF

	DO iz = 0,NZP-1
		DO ix = 0,NX/2
			UMY = UY(ix,0:NY-1,iz)
			CALL cfftf(NY,UMY,trigy)
			U(ix,0:NY/2-1,iz) = UMY(0:NY/2-1)/DFLOAT(NY)
			IF(INDEX_NXYZ==0) THEN
				U(ix,NY/2,iz) = DCMPLX( 0.0D0 , 0.0D0 )
			ELSE
				U(ix,NY/2,iz) = UMY(NY/2)/DFLOAT(NY+NY)
			ENDIF
			U(ix,NY/2+1:NY-1,iz) = UMY(NY/2+1:NY-1)/DFLOAT(NY)
		ENDDO
	ENDDO

	GOTO 30

! ---
20	DO iz = 0,NZP-1
		DO ix = 0,NX/2
			UMY = U(ix,0:NY-1,iz)
			UMY(NY/2) = 2.0D0*UMY(NY/2)
			CALL cfftb(NY,UMY,trigy)
			UY(ix,0:NY-1,iz) = UMY
		ENDDO
	ENDDO	  
	  
	IF(INDEX_Serial==1) THEN
		UK = UY
	ELSE
		DO K = 0,NPROC-1
			CALL MPI_SENDRECV(UY(0:NX/2,K*NYP:(K+1)*NYP-1,0:NZP-1),(NX/2+1)*NYP*NZP,MPI_COMPLEX16,K,10, &
							UK(0:NX/2,0:NYP-1,K*NZP:(K+1)*NZP-1),(NX/2+1)*NYP*NZP,MPI_COMPLEX16,K,10,MPI_COMM_WORLD,STATUS,ierr)
		ENDDO

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	ENDIF

	DO ix = 0,NX/2
		DO iy = 0,NYP-1
			UMZ = UK(ix,iy,0:NZ-1)
			UMZ(NZ/2) = 2.0D0*UMZ(NZ/2)
			CALL cfftb(NZ,UMZ,trigz)
			UK(ix,iy,0:NZ-1) = UMZ
		ENDDO
	ENDDO	  
	  
	DO iy = 0,NYP-1
		DO iz = 0,NZ-1
			UMX(0) = DREAL(UK(0,iy,iz))
			DO ix=1,NX/2-1
				UMX(2*ix-1) = DREAL(UK(ix,iy,iz))
				UMX(2*ix  ) = DIMAG(UK(ix,iy,iz))
			ENDDO
			UMX(NX-1) = 2.0D0*DREAL(UK(NX/2,iy,iz))
			CALL rfftb(NX,UMX,trigx)
			UP(0:NX-1,iy,iz) = UMX
		ENDDO
	ENDDO

! ---
30	DEALLOCATE(UMX,UMY,UMZ,UK,UY)

	RETURN
	END