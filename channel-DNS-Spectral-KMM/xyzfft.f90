! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE xyzfft (C,UP,U)

	  USE FFT

	  IMPLICIT REAL*8 (A-H, O-Z)

  include 'mpif.h'
  integer status(MPI_Status_Size)

      CHARACTER C
      REAL*8 UP(0:NX-1,0:NY,0:NZP-1)
	  COMPLEX*16 U(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)

	  REAL*8, ALLOCATABLE :: UMX(:),UMY(:),UK(:,:,:),UZ(:,:,:)
      COMPLEX*16, ALLOCATABLE :: UMZ(:)

      ALLOCATE(UMX(0:NX-1),UMY(0:NY),UMZ(0:NZ-1),UK(0:NX-1,0:NY,0:NZP-1),UZ(0:NXP-1,0:NY,0:NZ-1))

	  IF(C .EQ. 'F') GO TO 10
	  IF(C .EQ. 'B') GO TO 20

 10   DO iy = 0,NY
		 DO iz = 0,NZP-1
	        UMX = UP(0:NX-1,iy,iz)
		    CALL rfftf(NX,UMX,trigx)
	        UK(0,iy,iz) = UMX(0)/DFLOAT(NX)
	        UK(1,iy,iz) = 0.0D0
	        UK(2:NX-1,iy,iz) = UMX(1:NX-2)/DFLOAT(NX)
	     ENDDO
	  ENDDO

	  DO iz = 0,NZP-1
	     DO ix = 0,NX-1
	        UMY = UK(ix,0:NY,iz)
		    CALL cost (NY+1,UMY,trigy)
		    UMY(0) = UMY(0)/2.0D0
		    UMY(ny) = UMY(ny)/2.0D0
		    UK(ix,0:NY,iz) = UMY/DFLOAT(NY)
	     ENDDO
	  ENDDO	  

  DO K = 0,NPROC-1
	 id_send=MOD(my_id+k,NPROC)
	 id_recv=MOD(my_id+NPROC-k,NPROC)
	 CALL MPI_SendRecv(UK(id_send*NXP:(id_send+1)*NXP-1,0:NY,0:NZP-1),NXPNY1NZP,MPI_REAL8,id_send,10, &
		UZ(0:NXP-1,0:NY,id_recv*NZP:(id_recv+1)*NZP-1),NXPNY1NZP,MPI_REAL8,id_recv,10,MPI_COMM_WORLD,STATUS,ierr)
  ENDDO
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!  	  UZ=UK

	  DO ix = 0,NXP/2-1
		 DO iy = 0,NY
	        DO iz = 0,NZ-1
	           UMZ(iz) = DCMPLX(UZ(2*ix,iy,iz),UZ(2*ix+1,iy,iz))
			ENDDO
		    CALL cfftf(NZ,UMZ,trigz)
	        DO iz = -NZ/2+1,-1
			   U(ix,iy,iz) = UMZ(NZ+iz)/DFLOAT(NZ)
			ENDDO
	        DO iz = 0,NZ/2-1
			   U(ix,iy,iz) = UMZ(iz)/DFLOAT(NZ)
			ENDDO
	     ENDDO
	  ENDDO

	  GOTO 30


 20   DO ix = 0,NXP/2-1
	     DO iy = 0,NY
	        DO iz = 0,NZ/2-1
	           UMZ(iz) = U(ix,iy,iz)
			ENDDO
			UMZ(NZ/2) = DCMPLX(0.0D0,0.0D0)
	        DO iz = -NZ/2+1,-1
	           UMZ(NZ+iz) = U(ix,iy,iz)
			ENDDO
		    CALL cfftb(NZ,UMZ,trigz)
			DO iz = 0,NZ-1
	           UZ(2*ix,iy,iz) = DREAL(UMZ(iz))
	           UZ(2*ix+1,iy,iz) = DIMAG(UMZ(iz))
			ENDDO
	     ENDDO
	  ENDDO

  DO K = 0,NPROC-1
	 id_send=MOD(my_id+k,NPROC)
	 id_recv=MOD(my_id+NPROC-k,NPROC)
	 CALL MPI_SendRecv(UZ(0:NXP-1,0:NY,id_send*NZP:(id_send+1)*NZP-1),NXPNY1NZP,MPI_REAL8,id_send,10, &
		UK(id_recv*NXP:(id_recv+1)*NXP-1,0:NY,0:NZP-1),NXPNY1NZP,MPI_REAL8,id_recv,10,MPI_COMM_WORLD,STATUS,ierr)
  ENDDO
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!  	  UK=UZ

	  DO iz = 0,NZP-1
	     DO ix = 0,NX-1
	        UMY = UK(ix,0:NY,iz)
		    UMY(0) = 2.0D0*UMY(0)
		    UMY(ny) = 2.0D0*UMY(ny)
		    CALL cost (NY+1,UMY,trigy)
	 	    UK(ix,0:NY,iz) = UMY/2.0D0
	     ENDDO
	  ENDDO

      DO iy = 0,NY
         DO iz = 0,NZP-1
	        UMX(0) = UK(0,iy,iz)
	        UMX(1:NX-2) = UK(2:NX-1,iy,iz)
		    UMX(NX-1) = 0.0D0
		    CALL rfftb(NX,UMX,trigx)
	        UP(0:NX-1,iy,iz) = UMX
	     ENDDO
	  ENDDO

 30	  DEALLOCATE(UMX,UMY,UMZ,UK,UZ)

      RETURN
	  END
