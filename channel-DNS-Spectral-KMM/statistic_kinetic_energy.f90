
subroutine cal_k(u)
USE FFT
use RS_budget
INCLUDE 'dim.h'

!并行
include 'mpif.h'
integer status(MPI_Status_Size)
integer ierr
INTEGER COUNT            

complex*16,intent(in) :: u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)


REAL*8 PI
REAL*8 TEMP1,TEMP2,TEMP3
real*8,allocatable :: up(:,:,:,:)
allocate(up(3,0:NX-1,0:NY,0:NZP-1))

!把谱空间速度场变到物理空间
PI=DACOS(-1.0D0)
DO I=1,3
CALL xyzfft ('B',up(I,0:NX-1,0:NY,0:NZP-1),u(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
ENDDO


TEMP1=0.0D0
TEMP2=0.0D0
TEMP3=0.0D0

DO iy=0,NY
DO ix=0,NX-1
DO iz=0,NZP-1
TEMP1=TEMP1+up(1,ix,iy,iz)**2   !uu
TEMP2=TEMP2+up(2,ix,iy,iz)**2   !vv
TEMP3=TEMP3+up(3,ix,iy,iz)**2   !ww
enddo
enddo
enddo

! EK=EKP     !串行版本
!将所有进程的EK归纳到0进程中
COUNT=1
CALL MPI_REDUCE(TEMP1,EKT1,COUNT,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(TEMP2,EKT2,COUNT,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(TEMP3,EKT3,COUNT,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF(my_id.eq.0) then
	EKT1=EKT1/DFLOAT(NX*(NY+1)*NZ)
	EKT2=EKT2/DFLOAT(NX*(NY+1)*NZ)
	EKT3=EKT3/DFLOAT(NX*(NY+1)*NZ)
	WRITE(*,*)EKT1,EKT2,EKT3
ENDIF
deallocate(up)

end



