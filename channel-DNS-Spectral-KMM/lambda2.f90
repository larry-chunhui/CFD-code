subroutine lamda(u)
use FFT
INCLUDE 'dim.h'
complex*16,intent(in) :: u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)
character :: step*2
real*8,allocatable :: dud(:,:,:,:),QQ(:,:,:),up(:,:,:)
complex*16,allocatable :: temp(:,:,:,:)
real*8 w(3,3),s(3,3),ww(3,3),ss(3,3)

allocate(dud(9,0:NX-1,0:NY,0:NZP-1),QQ(0:NX-1,0:NY,0:NZP-1),up(0:NX-1,0:NY,0:NZP-1))
allocate(TEMP(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
write(STEP,'(i2)') my_id

CALL DERIV (1,TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !dud
CALL DERIV (2,TEMP(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i =1,9
CALL xyzfft ('B',dud(i,0:NX-1,0:NY,0:NZP-1),temp(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !dud
enddo

CALL xyzfft ('B',up(0:NX-1,0:NY,0:NZP-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))

DO iz=0,NZP-1
DO iy=0,NY
DO ix=0,NX-1
s(1,1)=0.5d0*(dud(1,ix,iy,iz)+dud(1,ix,iy,iz))
s(1,2)=0.5d0*(dud(4,ix,iy,iz)+dud(2,ix,iy,iz))
s(1,3)=0.5d0*(dud(7,ix,iy,iz)+dud(3,ix,iy,iz))
s(2,1)=s(1,2)
s(2,2)=0.5d0*(dud(5,ix,iy,iz)+dud(5,ix,iy,iz))
s(2,3)=0.5d0*(dud(8,ix,iy,iz)+dud(6,ix,iy,iz))
s(3,1)=s(1,3)
s(3,2)=s(2,3)
s(3,3)=0.5d0*(dud(9,ix,iy,iz)+dud(9,ix,iy,iz))
call mmultiple(ss,s,s)
w(1,1)=0.0d0
w(1,2)=0.5d0*(dud(2,ix,iy,iz)-dud(4,ix,iy,iz))
w(1,3)=0.5d0*(dud(3,ix,iy,iz)-dud(7,ix,iy,iz))
w(2,1)=-w(1,2)
w(2,2)=0.0d0
w(2,3)=0.5d0*(dud(6,ix,iy,iz)-dud(8,ix,iy,iz))
w(3,1)=-w(1,3)
w(3,2)=-w(2,3)
w(3,3)=0.0d0
call mmultiple(ww,w,w)
QQ(ix,iy,iz)=0.5d0*(dabs(ww(1,1))+dabs(ww(2,2))+dabs(ww(3,3))-dabs(ss(1,1))-dabs(ss(2,2))-dabs(ss(3,3)))
enddo
enddo
enddo

open(2,file='dud'//STEP//'.txt')
open(3,file='Q'//STEP//'.txt')
open(4,file='deltu'//STEP//'.txt')
do k=0,NZP-1
do j=0,NY
do i=0,NX-1
write(2,*) dud(1,i,j,k),dud(2,i,j,k),dud(3,i,j,k),dud(4,i,j,k),dud(5,i,j,k),dud(6,i,j,k),dud(7,i,j,k),dud(8,i,j,k),dud(9,i,j,k)
write(3,*) QQ(i,j,k)
! write(4,*) up(i,j,k)-aver1(1,i,j)
write(4,*) up(i,j,k)
enddo
enddo
enddo
close(2)
close(3)
close(4)

deallocate(dud,QQ)
deallocate(temp)
end

! subroutine lamda(u)
! use FFT
! ! use NUMERICAL_LIBRARIES
! INCLUDE 'dim.h'
! complex*16,intent(in) :: u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)
! character :: step*2
! real*8,allocatable :: dud(:,:,:,:)
! complex*16,allocatable :: temp(:,:,:,:)
! real*8 w(3,3),s(3,3),ww(3,3),ss(3,3),ssww(3,3)
! complex*16 lambda(3)
! allocate(dud(9,0:NX-1,0:NY,0:NZP-1))
! allocate(TEMP(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! write(STEP,'(i2)') my_id
! 
! CALL DERIV (1,TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !dud
! CALL DERIV (2,TEMP(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (3,TEMP(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (1,TEMP(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (2,TEMP(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (3,TEMP(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (1,TEMP(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (2,TEMP(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! CALL DERIV (3,TEMP(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! do i =1,9
! CALL xyzfft ('B',dud(i,0:NX-1,0:NY,0:NZP-1),temp(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !dud
! enddo
! 
! DO iz=0,NZP-1
! DO iy=0,NY
! DO ix=0,NX-1
! s(1,1)=0.5d0*(dud(1,ix,iy,iz)+dud(1,ix,iy,iz))
! s(1,2)=0.5d0*(dud(4,ix,iy,iz)+dud(2,ix,iy,iz))
! s(1,3)=0.5d0*(dud(7,ix,iy,iz)+dud(3,ix,iy,iz))
! s(2,1)=s(1,2)
! s(2,2)=0.5d0*(dud(5,ix,iy,iz)+dud(5,ix,iy,iz))
! s(2,3)=0.5d0*(dud(8,ix,iy,iz)+dud(6,ix,iy,iz))
! s(3,1)=s(1,3)
! s(3,2)=s(2,3)
! s(3,3)=0.5d0*(dud(9,ix,iy,iz)+dud(9,ix,iy,iz))
! call mmultiple(ss,s,s)
! w(1,1)=0.0d0
! w(1,2)=0.5d0*(dud(2,ix,iy,iz)-dud(4,ix,iy,iz))
! w(1,3)=0.5d0*(dud(3,ix,iy,iz)-dud(7,ix,iy,iz))
! w(2,1)=-w(1,2)
! w(2,2)=0.0d0
! w(2,3)=0.5d0*(dud(6,ix,iy,iz)-dud(8,ix,iy,iz))
! w(3,1)=-w(1,3)
! w(3,2)=-w(2,3)
! w(3,3)=0.0d0
! call mmultiple(ww,w,w)
! call msum(ssww,ss,ww)
! call devlrg(3,ssww,3,lambda)
! lambda2(ix,iy,iz)=DREAL(lambda(2))
! enddo
! enddo
! enddo
! 
! open(2,file='dud'//STEP//'.txt')
! do k=0,NZP-1
! do j=0,NY
! do i=0,NX-1
! write(2,*) 
! enddo
! enddo
! enddo
! close(2)
! 
! deallocate(dud)
! deallocate(temp)
! end
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine mmultiple(MAT,MAT1,MAT2)
      implicit none
	  integer i,j
      real*8 MAT(3,3),MAT1(3,3),MAT2(3,3)

      do j=1,3
      do i=1,3
      MAT(i,j)=MAT1(i,1)*MAT2(1,j)+MAT1(i,2)*MAT2(2,j)+MAT1(i,3)*MAT2(3,j)
      enddo
      enddo
      
      return
      end subroutine


      subroutine msum(MAT,MAT1,MAT2)
      implicit none
	  integer i,j      
      real*8 MAT(3,3),MAT1(3,3),MAT2(3,3)

      do j=1,3
      do i=1,3
      MAT(i,j)=MAT1(i,j)+MAT2(i,j)
      enddo
      enddo
      
      return
      end subroutine
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@