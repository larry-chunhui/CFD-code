
subroutine plas_actuator
use plas_v
use FFT
include 'dim.h'
! include 'plas_variable.h'
integer :: I0,m
real*8 :: x1,y1,tempE,kx,ky
REAL*8,ALLOCATABLE :: plas_E(:,:,:),plas_Fp(:,:,:,:)
REAL*8,ALLOCATABLE :: x(:,:,:),y(:,:,:),z(:,:,:)
real*8 :: temp_Fp(2,0:NX-1,0:NY,0:NZP-1),temp1_Fp(2,0:NX-1,0:NY,0:NZP-1),temp2_Fp(2,0:NX-1,0:NY,0:NZP-1)
real*8 :: PI,Dc,ref_E,ref_L,plas_time,plas_Eb,plas_H,plas_W,plas_S,plas_fq,plas_dt,kksi,keta
common /plas_position/I0
common /ref/ref_E,ref_L
common /acttime/plas_time
real*8 :: plas_a,plas_b,plas_d,plas_k1,plas_k2,plas_k3
common /plas_para/Dc,plas_a,plas_b,plas_d,plas_fq,plas_dt,plas_k1,plas_k2,kksi,keta
CHARACTER STEP*2
write(STEP,'(i2)') my_id

ALLOCATE(x(0:NX-1,0:NY,0:NZP-1),y(0:NX-1,0:NY,0:NZP-1),z(0:NX-1,0:NY,0:NZP-1))
ALLOCATE(plas_E(0:NX-1,0:NY,0:NZP-1),plas_Fp(2,0:NX-1,0:NY,0:NZP-1))


PI = 2.0D0*DACOS(0.0D0)
do k=0,NZP-1
do j=0,NY
do i=0,NX-1
   x(i,j,k)=8*PI*DFLOAT(i)/DFLOAT(NX)
   y(i,j,k)=DCOS(PI*DFLOAT(j)/DFLOAT(NY))
   z(i,j,k)=2*PI*DFLOAT(my_id*NZP+k)/DFLOAT(NZ)
enddo
enddo
enddo

open(1,file='plasma.txt',status='unknown',action='read')
read(1,*) I0
read(1,*) Dc
read(1,*) ref_E
read(1,*) ref_L
read(1,*) plas_time
read(1,*) plas_Eb
read(1,*) plas_fq
read(1,*) plas_dt
read(1,*) plas_H
read(1,*) plas_W
read(1,*) plas_S
read(1,*) kksi
read(1,*) keta
close(1)

write(*,*) 'plasma parameter:'
write(*,*) 'actuator position =',x(I0,0,0)
write(*,*) 'plasma force scale =',Dc
write(*,*) 'discharge time =',plas_dt
write(*,*) 'breakdown field strength =',plas_Eb
write(*,*) 'plasma separation =',plas_S,'height =',plas_H,'width =',plas_W

!non_dimension
plas_Eb=plas_Eb/ref_E
plas_a=plas_H/ref_L
plas_b=plas_W/ref_L
plas_d=plas_S/ref_L
plas_k1=(1.0d0-plas_Eb)/plas_b
plas_k2=(1.0d0-plas_Eb)/plas_a
kx=1/dsqrt(1+(plas_k1/plas_k2)**2.0d0)
ky=1/dsqrt(1+(plas_k2/plas_k1)**2.0d0)

plas_k3=(1.0d0-plas_Eb)/plas_d

!electric field and plas_force
do k=0,NZP-1
do j=0,NY/2-1
do i=0,NX-1
   x1=x(i,j,k)-x(I0,j,k)
   y1=1-y(i,j,k)
   plas_E(i,j,k)=1-plas_k1*x1-plas_k2*y1

   if (x1 .lt. 0.0d0)  plas_E(i,j,k)=1-plas_k3*dabs(x1)-plas_k2*y1

   tempE=plas_E(i,j,k)
!    if(x1 .lt. 0.0D0 .OR. x1 .gt. plas_b .OR. y1 .gt. plas_a) plas_E(i,j,k)=0.0d0
   if(x1 .gt. plas_b .OR. y1 .gt. plas_a) plas_E(i,j,k)=0.0d0
   if(tempE .lt. plas_Eb) plas_E(i,j,k)=0.0d0
   plas_Fp(1,i,j,k)=plas_fq*plas_dt*Dc*kx*plas_E(i,j,k)*kksi
   plas_Fp(2,i,j,k)=plas_fq*plas_dt*Dc*ky*plas_E(i,j,k)*keta
   plas_Fp(1,i,NY-j,k)=plas_Fp(1,i,j,k)
   plas_Fp(2,i,NY-j,k)=-plas_Fp(2,i,j,k)
enddo
enddo
enddo
plas_Fp(:,:,NY/2,:)=0.0d0

DO I=1,2
   CALL xyzfft ('F',plas_Fp(I,0:NX-1,0:NY,0:NZP-1),plas_F(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
ENDDO
plas_F(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)=0.0d0

! output body force
open(1,file='plas_force'//STEP//'.plt')
WRITE(1,130)'"X" ,','"Y" ,','"Z" ,','"Fpx" ,','"Fpy" ',NX,NY+1,NZP
do k=0,NZP-1
do j=0,NY
do i=0,NX-1
   write(1,*) x(i,j,k),y(i,j,k),z(i,j,k),plas_Fp(1,i,j,k),plas_Fp(2,i,j,k)
enddo
enddo
enddo
close(1)
130  FORMAT('VARIABLES = ',5(1X,A)/'ZONE ,  I=',I3,',  J=',I3,', k=',I3,',  F=POINT')

END
!####################################################
! do m=1,2
! do k=0,NZP-1
! do j=0,NY
! do i=0,NX-1
! temp2_Fp(m,i,j,k)=plas_Fp(m,i ,0,k)
! enddo
! enddo
! enddo
! enddo
! 
! open(1,file='temp2_Fp_sy'//STEP//'.plt')
! WRITE(1,132)'"X" ,','"Y" ,','"Fpx" ,','"Fpy" ',NX,NY+1
! do k=4,4
! do j=0,NY
! do i=0,NX-1
! write(1,*) x(i,j,k),y(i,j,k),temp2_Fp(1,i,j,k),temp2_Fp(2,i,j,k)
! enddo
! enddo
! enddo
! close(1)
! 132  FORMAT('VARIABLES = ',4(1X,A)/'ZONE ,  I=',I3,',  J=',I3,',  F=POINT')
! 
! !check FFT direction X temp2_Fp
! DO I=1,2
!    CALL xyzfft ('F',temp2_Fp(I,0:NX-1,0:NY,0:NZP-1),plas_F(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! ENDDO
! 
! DO I=1,2
!    CALL xyzfft ('B',temp_Fp(I,0:NX-1,0:NY,0:NZP-1),plas_F(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! ENDDO
! 
! temp2_Fp=temp2_Fp-temp_Fp
! 
! open(1,file='temp2_err_sy'//STEP//'.plt')
! WRITE(1,131)'"X" ,','"Y" ,','"Fpx" ,','"Fpy" ',NX,NY+1
! do k=4,4
! do j=0,NY
! do i=0,NX-1
!    write(1,*) x(i,j,k),y(i,j,k),temp2_Fp(1,i,j,k),temp2_Fp(2,i,j,k)
! enddo
! enddo
! enddo
! close(1)
! 131  FORMAT('VARIABLES = ',4(1X,A)/'ZONE ,  I=',I3,',  J=',I3,',  F=POINT')

