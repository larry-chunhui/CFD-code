! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE xyzfft2 (C,UP2,U)

	  USE FFT

	  IMPLICIT REAL*8 (A-H, O-Z)

  include 'mpif.h'
  integer status(MPI_Status_Size)

      CHARACTER C
      REAL*8 UP2(0:NX2-1,0:NY2,0:NZ2P-1)
	  COMPLEX*16 U(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)

	  REAL*8, ALLOCATABLE :: UMX(:),UMY(:),UK(:,:,:),UZ(:,:,:)
      COMPLEX*16, ALLOCATABLE :: UMZ(:)

      ALLOCATE(UMX(0:NX2-1),UMY(0:NY2),UMZ(0:NZ2-1),UK(0:NX-1,0:NY2,0:NZ2P-1),UZ(0:NXP-1,0:NY,0:NZ2-1))

	  IF(C .EQ. 'F') GO TO 10
	  IF(C .EQ. 'B') GO TO 20

 10   DO iy = 0,NY2
         DO iz = 0,NZ2P-1
	        UMX = UP2(0:NX2-1,iy,iz)
		    CALL rfftf(NX2,UMX,trigx2)
	        UK(0,iy,iz) = UMX(0)/DFLOAT(NX2)
	        UK(1,iy,iz) = 0.0D0
	        UK(2:NX-1,iy,iz) = UMX(1:NX-2)/DFLOAT(NX2)
	     ENDDO
	  ENDDO

	  DO iz = 0,NZ2P-1
	     DO ix = 0,NX-1
	        UMY = UK(ix,0:NY2,iz)
		    CALL cost (NY2+1,UMY,trigy2)
		    UMY(0) = UMY(0)/2.0D0
		    UMY(NY2) = UMY(NY2)/2.0D0
		    UK(ix,0:NY,iz) = UMY(0:NY)/DFLOAT(NY2)
	     ENDDO
	  ENDDO

  DO K = 0,NPROC-1
	 id_send=MOD(my_id+k,NPROC)
	 id_recv=MOD(my_id+NPROC-k,NPROC)
	 CALL MPI_SendRecv(UK(id_send*NXP:(id_send+1)*NXP-1,0:NY,0:NZ2P-1),NXPNY1NZ2P,MPI_REAL8,id_send,10, &
		UZ(0:NXP-1,0:NY,id_recv*NZ2P:(id_recv+1)*NZ2P-1),NXPNY1NZ2P,MPI_REAL8,id_recv,10,MPI_COMM_WORLD,STATUS,ierr)
  ENDDO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  	  UZ=UK(0:NXP-1,0:NY,0:NZ2-1)

	  DO ix = 0,NXP/2-1
	     DO iy = 0,NY
	        DO iz = 0,NZ2-1
	           UMZ(iz) = DCMPLX(UZ(2*ix,iy,iz),UZ(2*ix+1,iy,iz))
			ENDDO
		    CALL cfftf(NZ2,UMZ,trigz2)
	        DO iz = -NZ/2+1,-1
			   U(ix,iy,iz) = UMZ(NZ2+iz)/DFLOAT(NZ2)
			ENDDO
	        DO iz = 0,NZ/2-1
			   U(ix,iy,iz) = UMZ(iz)/DFLOAT(NZ2)
			ENDDO
	     ENDDO
	  ENDDO

	  GOTO 30


 20	  DO ix = 0,NXP/2-1
	     DO iy = 0,NY
	        DO iz = 0,NZ/2-1
	           UMZ(iz) = U(ix,iy,iz)
			ENDDO
	        UMZ(NZ/2:NZ2-NZ/2) = DCMPLX(0.0D0,0.0D0)
	        DO iz = -NZ/2+1,-1
	           UMZ(NZ2+iz) = U(ix,iy,iz)
			ENDDO
		    CALL cfftb(NZ2,UMZ,trigz2)
			DO iz = 0,NZ2-1
	           UZ(2*ix,iy,iz) = DREAL(UMZ(iz))
	           UZ(2*ix+1,iy,iz) = DIMAG(UMZ(iz))
			ENDDO
	     ENDDO
	  ENDDO

  DO K = 0,NPROC-1
	 id_send=MOD(my_id+k,NPROC)
	 id_recv=MOD(my_id+NPROC-k,NPROC)

	 CALL MPI_SendRecv(UZ(0:NXP-1,0:NY,id_send*NZ2P:(id_send+1)*NZ2P-1),NXPNY1NZ2P,MPI_REAL8,id_send,10, &
		UK(id_recv*NXP:(id_recv+1)*NXP-1,0:NY,0:NZ2P-1),NXPNY1NZ2P,MPI_REAL8,id_recv,10,MPI_COMM_WORLD,STATUS,ierr)
  ENDDO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!  	  UK(0:NX-1,0:NY,0:NZ2P-1)=UZ

	  DO iz = 0,NZ2P-1
	     DO ix = 0,NX-1
	        UMY(NY+1:NY2) = 0.0D0
	        UMY(0:NY) = UK(ix,0:NY,iz)
		    UMY(0) = 2.0D0*UMY(0)
		    UMY(NY2) = 2.0D0*UMY(NY2)
		    CALL cost (NY2+1,UMY,trigy2)
	 	    UK(ix,0:NY2,iz) = UMY/2.0D0
	     ENDDO
	  ENDDO

      DO iy = 0,NY2
         DO iz = 0,NZ2P-1
	        UMX(0) = UK(0,iy,iz)
	        UMX(1:NX-2) = UK(2:NX-1,iy,iz)
		    UMX(NX-1:NX2-1) = 0.0D0
		    CALL rfftb(NX2,UMX,trigx2)
	        UP2(0:NX2-1,iy,iz) = UMX
	     ENDDO
	  ENDDO

 30	  DEALLOCATE(UMX,UMY,UMZ,UK,UZ)

      RETURN
	  END
