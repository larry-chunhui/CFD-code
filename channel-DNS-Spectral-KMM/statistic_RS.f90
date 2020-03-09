
subroutine RS_sum(u,p)
USE FFT
use RS_budget
INCLUDE 'dim.h'
complex*16,intent(in) :: u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),p(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)
integer :: il,ik,im
REAL*8 PI
real*8,allocatable :: up(:,:,:,:),pp(:,:,:),uu(:,:,:,:),uuu(:,:,:,:),uuuu(:,:,:,:),&
					  dpd(:,:,:,:),dud(:,:,:,:),duud(:,:,:,:),&
                      duuud(:,:,:,:),ddudd(:,:,:,:),dduudd(:,:,:,:),wp(:,:,:,:)
real*8,allocatable :: dpd1(:,:,:,:)
complex*16,allocatable :: temp(:,:,:),temp1(:,:,:,:),temp2(:,:,:,:),temp3(:,:,:,:)
allocate(up(3,0:NX-1,0:NY,0:NZP-1),wp(3,0:NX-1,0:NY,0:NZP-1),pp(0:NX-1,0:NY,0:NZP-1),&
		 uu(6,0:NX-1,0:NY,0:NZP-1),dpd(3,0:NX-1,0:NY,0:NZP-1),&
         dud(9,0:NX-1,0:NY,0:NZP-1),duud(16,0:NX-1,0:NY,0:NZP-1),duuud(12,0:NX-1,0:NY,0:NZP-1),&
         uuu(3,0:NX-1,0:NY,0:NZP-1),uuuu(3,0:NX-1,0:NY,0:NZP-1),ddudd(9,0:NX-1,0:NY,0:NZP-1),dduudd(12,0:NX-1,0:NY,0:NZP-1))
allocate(dpd1(5,0:NX-1,0:NY,0:NZP-1))
allocate(temp(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp1(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),&
		 temp2(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp3(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
!把谱空间速度、压力场变到物理空间
PI=DACOS(-1.0D0)
DO I=1,3
CALL xyzfft ('B',up(I,0:NX-1,0:NY,0:NZP-1),u(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
ENDDO
CALL xyzfft ('B',pp(0:NX-1,0:NY,0:NZP-1),p(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
DO iy=0,NY
DO ix=0,NX-1
DO iz=0,NZP-1
uu(1,ix,iy,iz)=up(1,ix,iy,iz)**2   !uu
uu(2,ix,iy,iz)=up(2,ix,iy,iz)**2   !vv
uu(3,ix,iy,iz)=up(3,ix,iy,iz)**2   !ww
uu(4,ix,iy,iz)=up(1,ix,iy,iz)*up(2,ix,iy,iz)   !uv
uu(5,ix,iy,iz)=up(1,ix,iy,iz)*up(3,ix,iy,iz)   !uw
uu(6,ix,iy,iz)=up(2,ix,iy,iz)*up(3,ix,iy,iz)   !vw

uuu(1,ix,iy,iz)=up(1,ix,iy,iz)**3   !uuu
uuu(2,ix,iy,iz)=up(2,ix,iy,iz)**3   !vvv
uuu(3,ix,iy,iz)=up(3,ix,iy,iz)**3   !www
uuuu(1,ix,iy,iz)=up(1,ix,iy,iz)**4  !uuuu
uuuu(2,ix,iy,iz)=up(2,ix,iy,iz)**4  !vvvv
uuuu(3,ix,iy,iz)=up(3,ix,iy,iz)**4  !wwww
pp(ix,iy,iz)=pp(ix,iy,iz)-(uu(1,ix,iy,iz)+uu(2,ix,iy,iz)+uu(3,ix,iy,iz))/2.0d0
enddo
enddo
enddo

CALL xyzfft ('F',pp(0:NX-1,0:NY,0:NZP-1),TEMP(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i=1,3
call DERIV (i,TEMP1(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !dpd
CALL xyzfft ('B',dpd(i,0:NX-1,0:NY,0:NZP-1),TEMP1(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))   !dpd
enddo

CALL DERIV (1,TEMP2(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !dud
CALL DERIV (2,TEMP2(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP2(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP2(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP2(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP2(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP2(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP2(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP2(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),u(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i =1,9
CALL xyzfft ('B',dud(i,0:NX-1,0:NY,0:NZP-1),temp2(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !dud
enddo
!计算涡量 10.24
DO iy=0,NY
DO ix=0,NX-1
DO iz=0,NZP-1
	wp(1,ix,iy,iz)=dud(8,ix,iy,iz)-dud(6,ix,iy,iz)
	wp(2,ix,iy,iz)=dud(3,ix,iy,iz)-dud(7,ix,iy,iz)
	wp(3,ix,iy,iz)=dud(4,ix,iy,iz)-dud(2,ix,iy,iz)
ENDDO
ENDDO
ENDDO
CALL DERIV (1,TEMP3(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !ddudd
CALL DERIV (2,TEMP3(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP3(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP3(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP3(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP3(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP3(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP3(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP3(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP2(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i =1,9
CALL xyzfft ('B',ddudd(i,0:NX-1,0:NY,0:NZP-1),temp3(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !ddudd
enddo

call dxyz(up,pp,uu,dpd1,duud,duuud,dduudd)  !dpdx,dpudx,duudx,duuudx,dduudxx

! sum
DO ix=0,NX-1
DO iz=0,NZP-1
sum0(1)=sum0(1)+up(1,ix,NY/2,iz)           !Uc
ENDDO
ENDDO

DO iy=0,NY
DO ix=0,NX-1
DO iz=0,NZP-1
sum0(2)=sum0(2)+up(1,ix,iy,iz)*(dcos(dfloat(iy)*PI/dfloat(NY))-&
                                dcos(dfloat(iy+1)*PI/dfloat(NY)))         !Um
sum1(1,iy)=sum1(1,iy)+up(1,ix,iy,iz)    !U
sum1(2,iy)=sum1(2,iy)+up(2,ix,iy,iz)    !V
sum1(3,iy)=sum1(3,iy)+up(3,ix,iy,iz)    !W

sum2(1,iy)=sum2(1,iy)+uu(1,ix,iy,iz)   !UU
sum2(2,iy)=sum2(2,iy)+uu(2,ix,iy,iz)   !VV
sum2(3,iy)=sum2(3,iy)+uu(3,ix,iy,iz)   !ww
sum2(4,iy)=sum2(4,iy)+uu(4,ix,iy,iz)   !uv
sum2(5,iy)=sum2(5,iy)+uu(5,ix,iy,iz)   !uw
sum2(6,iy)=sum2(6,iy)+uu(6,ix,iy,iz)   !vw

sum3(1,iy)=sum3(1,iy)+uuu(1,ix,iy,iz)  !uuu
sum3(2,iy)=sum3(2,iy)+uuu(2,ix,iy,iz)  !vvv
sum3(3,iy)=sum3(3,iy)+uuu(3,ix,iy,iz)  !www

sum4(1,iy)=sum4(1,iy)+uuuu(1,ix,iy,iz)  !uuuu
sum4(2,iy)=sum4(2,iy)+uuuu(2,ix,iy,iz)  !vvvv
sum4(3,iy)=sum4(3,iy)+uuuu(3,ix,iy,iz)  !wwww

sumwp(1,iy)=sumwp(1,iy)+wp(1,ix,iy,iz)  !Wx
sumwp(2,iy)=sumwp(2,iy)+wp(2,ix,iy,iz)  !Wy
sumwp(3,iy)=sumwp(3,iy)+wp(3,ix,iy,iz)  !Wz

sumwp2(1,iy)=sumwp2(1,iy)+wp(1,ix,iy,iz)**2  !Wx*Wx
sumwp2(2,iy)=sumwp2(2,iy)+wp(2,ix,iy,iz)**2  !Wy*Wy
sumwp2(3,iy)=sumwp2(3,iy)+wp(3,ix,iy,iz)**2  !Wz*Wz

sump(1,iy)=sump(1,iy)+pp(ix,iy,iz)    !p
sump(2,iy)=sump(2,iy)+pp(ix,iy,iz)**2.0d0   !p*p
sump(3,iy)=sump(3,iy)+dpd(1,ix,iy,iz)       !dpdx
sump(4,iy)=sump(4,iy)+dpd(2,ix,iy,iz)       !dpdy
sump(5,iy)=sump(5,iy)+dpd(3,ix,iy,iz)       !dpdz
! sump(6,ix,iy)=sump(6,ix,iy)+dpd1(1,ix,iy,iz)       !dpudx
! sump(7,ix,iy)=sump(7,ix,iy)+dpd1(2,ix,iy,iz)       !dpudy
! sump(8,ix,iy)=sump(8,ix,iy)+dpd1(3,ix,iy,iz)       !dpvdx
! sump(9,ix,iy)=sump(9,ix,iy)+dpd1(4,ix,iy,iz)       !dpvdy
! sump(10,ix,iy)=sump(10,ix,iy)+dpd1(5,ix,iy,iz)     !dpwdz
! sump(11,ix,iy)=sump(11,ix,iy)+pp(ix,iy,iz)*dud(1,ix,iy,iz) !p*dudx
! sump(12,ix,iy)=sump(12,ix,iy)+pp(ix,iy,iz)*dud(2,ix,iy,iz) !p*dudy
! sump(13,ix,iy)=sump(13,ix,iy)+pp(ix,iy,iz)*dud(4,ix,iy,iz) !p*dvdx
! sump(14,ix,iy)=sump(14,ix,iy)+pp(ix,iy,iz)*dud(5,ix,iy,iz) !p*dvdy
! sump(15,ix,iy)=sump(15,ix,iy)+pp(ix,iy,iz)*dud(9,ix,iy,iz) !p*dwdz
sump(6,iy)=sump(6,iy)+dpd(1,ix,iy,iz)*up(1,ix,iy,iz)  !u*dpdx
sump(7,iy)=sump(7,iy)+dpd(2,ix,iy,iz)*up(2,ix,iy,iz)  !v*dpdy
sump(8,iy)=sump(8,iy)+dpd(3,ix,iy,iz)*up(3,ix,iy,iz)  !w*dpdz
sump(9,iy)=sump(9,iy)+dpd(2,ix,iy,iz)*up(1,ix,iy,iz)  !u*dpdy
sump(10,iy)=sump(10,iy)+dpd(1,ix,iy,iz)*up(2,ix,iy,iz)  !v*dpdx

do il=1,16
sumduud(il,iy)=sumduud(il,iy)+duud(il,ix,iy,iz)  !duud
enddo

do ik=1,9
sumdud(ik,iy)=sumdud(ik,iy)+dud(ik,ix,iy,iz)    !dud
sumdd(ik,iy)=sumdd(ik,iy)+dud(ik,ix,iy,iz)**2.0d0  !dudx*dudx
sumddudd(ik,iy)=sumddudd(ik,iy)+ddudd(ik,ix,iy,iz)  !ddudd
enddo

sumdd(10,iy)=sumdd(10,iy)+dud(1,ix,iy,iz)*dud(4,ix,iy,iz)  !dudx*dvdx
sumdd(11,iy)=sumdd(11,iy)+dud(2,ix,iy,iz)*dud(5,ix,iy,iz)  !dudy*dvdy
sumdd(12,iy)=sumdd(12,iy)+dud(3,ix,iy,iz)*dud(6,ix,iy,iz)  !dudz*dvdz

do im=1,12
sumdduudd(im,iy)=sumdduudd(im,iy)+dduudd(im,ix,iy,iz) !dduudd
sumduuud(im,iy)=sumduuud(im,iy)+duuud(im,ix,iy,iz) !duuud
enddo

enddo
enddo
enddo

deallocate(up,pp,uu,dpd,dud,duud,duuud,ddudd,dduudd)
deallocate(dpd1)
deallocate(temp,temp1,temp2,temp3)
end


subroutine RS_statis(NNN)
USE FFT
use RS_budget
INCLUDE 'dim.h'
integer,intent(in) :: NNN
integer :: I0
real*8 :: CF_sum
character :: step*5
character :: step2*2
integer :: STI(4)
real*8  :: tao0,utao
real*8,allocatable :: xx(:),yy(:),yp(:),prms(:),averu(:),shearF(:),uu(:,:)
real*8,allocatable :: urms(:),vrms(:),wrms(:),uv(:),wprms(:,:)
real*8 PI
real*8,allocatable ::S(:,:),F(:,:)
!*************************
!暂时先忽略K方程
!*************************
! REAL*8,ALLOCATABLE :: P_11(:,:),T_11(:,:),D_11(:,:),PS_11(:,:),PD_11(:,:),ipsron_11(:,:),&
!                       P_22(:,:),T_22(:,:),D_22(:,:),PS_22(:,:),PD_22(:,:),ipsron_22(:,:),&
!                       P_33(:,:),T_33(:,:),D_33(:,:),PS_33(:,:),PD_33(:,:),ipsron_33(:,:),&
! 					  P_12(:,:),T_12(:,:),D_12(:,:),PS_12(:,:),PD_12(:,:),ipsron_12(:,:),&
! 					  P_k(:,:) ,T_k(:,:) ,D_k(:,:) ,PD_k(:,:) ,ipsron_k(:,:)
! ALLOCATE(P_11(0:NX-1,0:NY),T_11(0:NX-1,0:NY),D_11(0:NX-1,0:NY),PS_11(0:NX-1,0:NY),&
!          PD_11(0:NX-1,0:NY),ipsron_11(0:NX-1,0:NY))
! ALLOCATE(P_22(0:NX-1,0:NY),T_22(0:NX-1,0:NY),D_22(0:NX-1,0:NY),PS_22(0:NX-1,0:NY),&
!          PD_22(0:NX-1,0:NY),ipsron_22(0:NX-1,0:NY))
! ALLOCATE(P_33(0:NX-1,0:NY),T_33(0:NX-1,0:NY),D_33(0:NX-1,0:NY),PS_33(0:NX-1,0:NY),&
!          PD_33(0:NX-1,0:NY),ipsron_33(0:NX-1,0:NY))
! ALLOCATE(P_12(0:NX-1,0:NY),T_12(0:NX-1,0:NY),D_12(0:NX-1,0:NY),PS_12(0:NX-1,0:NY),&
!          PD_12(0:NX-1,0:NY),ipsron_12(0:NX-1,0:NY))
! ALLOCATE(P_k(0:NX-1,0:NY),T_k(0:NX-1,0:NY),D_k(0:NX-1,0:NY),PD_k(0:NX-1,0:NY),ipsron_k(0:NX-1,0:NY))

allocate(xx(0:NX-1),yy(0:NY),yp(0:NY/2),prms(0:NY),&
         averu(0:NY),shearF(0:NY),uu(6,0:NY))
ALLOCATE(urms(0:NY),vrms(0:NY),wrms(0:NY),uv(0:NY),wprms(3,0:NY))
ALLOCATE(S(1:3,0:NY),F(1:3,0:NY))

STI(1)=NX/8.0d0;STI(2)=3.0d0*NX/8.0d0;STI(3)=5.0d0*NX/8.0d0;STI(4)=7.0d0*NX/8.0d0


write(STEP,'(i5)') NNN


PI = 2.0D0*DACOS(0.0D0)
I0=64


tao0=(dabs(averdud(2,0))+dabs(averdud(2,NY)))/2.0d0
utao=dsqrt(tao0/Re)
! utao(ix)=dsqrt(tao0(ix))      CHANGE BY LIU CHUN HUI

do ix=0,NX-1
xx(ix)=4*PI*DFLOAT(ix)/DFLOAT(NX)
enddo
do iy=0,NY
yy(iy)=DCOS(PI*DFLOAT(iy)/DFLOAT(NY))
enddo

DO iy=0,NY/2
! yp(ix,iy)=(1.0d0-yy(iy))*dsqrt(Re)*utao(ix)       CHANGE BY LIU CHUN HUI
yp(iy)=(1.0d0-yy(iy))*Re*utao
enddo

DO iy=0,NY
uu(1,iy)=aver2(1,iy)-aver1(1,iy)**2.0d0  !uu
uu(2,iy)=aver2(2,iy)-aver1(2,iy)**2.0d0  !vv
uu(3,iy)=aver2(3,iy)-aver1(3,iy)**2.0d0  !ww
uu(4,iy)=aver2(4,iy)-aver1(1,iy)*aver1(2,iy)  !uv
uu(5,iy)=aver2(5,iy)-aver1(1,iy)*aver1(3,iy)  !uw
uu(6,iy)=aver2(6,iy)-aver1(2,iy)*aver1(3,iy)  !vw
enddo

DO iy=0,NY
! averu(ix,iy)=aver1(1,ix,iy)*dsqrt(1.0d0/tao0(ix))*dsqrt(Re)
! urms(ix,iy)=dsqrt(dabs(uu(1,ix,iy)))
! vrms(ix,iy)=dsqrt(dabs(uu(2,ix,iy)))
! wrms(ix,iy)=dsqrt(dabs(uu(3,ix,iy)))
! uv(ix,iy)=uu(4,ix,iy)
! prms(ix,iy)=averp(2,ix,iy)-averp(1,ix,iy)**2.0d0
! prms(ix,iy)=dsqrt(dabs(prms(ix,iy)))
! shearF(ix,iy)=averdud(2,ix,iy)/Re-uu(4,ix,iy)

!gui yi hua by ut
averu(iy)=aver1(1,iy)/utao
urms(iy)=dsqrt(dabs(uu(1,iy)))/utao
vrms(iy)=dsqrt(dabs(uu(2,iy)))/utao
wrms(iy)=dsqrt(dabs(uu(3,iy)))/utao

wprms(1,iy)=dsqrt(averwp2(1,iy)-averwp(1,iy)**2)/(Re*utao**2)  !Wxrms
wprms(2,iy)=dsqrt(averwp2(2,iy)-averwp(2,iy)**2)/(Re*utao**2)  !Wyrms
wprms(3,iy)=dsqrt(averwp2(3,iy)-averwp(3,iy)**2)/(Re*utao**2)  !Wzrms

uv(iy)=uu(4,iy)/tao0
prms(iy)=averp(2,iy)-averp(1,iy)**2.0d0
prms(iy)=dsqrt(dabs(prms(iy)))/(utao**2)
shearF(iy)=(averdud(2,iy)/Re-uu(4,iy))/(utao**2)
!ADD BY LIU CHUN HUI
!Compute flateness and Skewness
!uu(1,iy)=aver2(1,iy)-aver1(1,iy)**2.0d0  !<u'u'>=<uu>-<u>*<u>
do j=1,3
	S(j,iy)=aver3(j,iy)-3.0d0*aver2(j,iy)*aver1(j,iy)+2.0d0*(aver1(j,iy)**3)
	S(j,iy)=S(j,iy)/(uu(j,iy)**(1.5d0))
	F(j,iy)=aver4(j,iy)-4.0d0*aver3(j,iy)*aver1(j,iy)+6.0d0*aver2(j,iy)*aver1(j,iy)**2-3.0d0*aver1(j,iy)**4
	F(j,iy)=F(j,iy)/(uu(j,iy)**2)
enddo

enddo


!************************
!BY LIU CHUNHUI 
!************************
!CX XU PIC.3.6 
OPEN(1,FILE='loglaw up.plt')
WRITE(1,*)'VARIABLES="y+","u+","Urms","Vrms","Wrms"'
DO iy=0,NY/2
	WRITE(1,*) yp(iy),averu(iy),urms(iy),vrms(iy),wrms(iy)
ENDDO
CLOSE(1)

!CX XU PIC.3.7A
OPEN(2,FILE='UVWrms Profile.plt')
WRITE(2,*)'VARIABLES="y","Urms","Vrms","Wrms","Prms","averU",",shearF","1/Redu/dy","<uv>","Wxrms","Wyrms","Wzrms"'
DO iy=0,NY
	WRITE(2,*) yy(iy),urms(iy),vrms(iy),wrms(iy),prms(iy),aver1(1,iy),shearF(iy),averdud(2,iy)/(Re*utao**2),&
	           -uu(4,iy)/(utao**2),wprms(1,iy),wprms(2,iy),wprms(3,iy)
ENDDO
CLOSE(2)

!KIM MOIN TABLE1
OPEN(3,FILE='FLOW VARIBLES.txt')
WRITE(3,*)'Uc=',aver0(1)
WRITE(3,*)'Um=',aver0(2)
WRITE(3,*)'Utao=',utao
CLOSE(3)

!CX XU PIC.3.11
OPEN(4,FILE='Flatness_Skewness.plt')
WRITE(4,*)'VARIABLES="y","Su","Sv","Sw","Fu","Fv","Fw"'
DO iy=1,NY-1
	WRITE(4,*)yy(iy),S(1,iy),S(2,iy),S(3,iy),F(1,iy),F(2,iy),F(3,iy)
ENDDO
CLOSE(4)




!sum1(1,iy)=sum1(1,iy)+up(1,ix,iy,iz)    !U
!************************
!BY HUANG JUN JI FOUR STATION PLOT
!************************

! open(1,file='aver_flowfield'//STEP//'.plt')
! write(1,130) '"x" ,','"y" ,','"p" ,','"U" ,','"V" ,','"W" ,','"UU" ,','"VV" ,','"WW" ,','"UV" ,','"UW" ,','"VW" ',NX,NY+1
! do iy=0,NY
! do ix=0,NX-1
! write(1,*) xx(ix),yy(iy),averp(1,ix,iy),aver1(1,ix,iy),aver1(2,ix,iy),aver1(3,ix,iy),&
!            -uu(1,ix,iy),-uu(2,ix,iy),-uu(3,ix,iy),-uu(4,ix,iy),-uu(5,ix,iy),-uu(6,ix,iy)
! enddo
! enddo
! 130 format('variables = ',12(1x,A)/'zone, I =',I3,', J =',I3,', F=POINT')
! close(1)
! 
! CF_sum=0.0d0
! open(2,file='tao0-x'//step//'.plt')
! write(2,142) '"x" ,','"Cf" ,','"averdp" ,',NX
! do ix=0,NX-1
! write(2,*) xx(ix),2*tao0(ix)/Re,averp(3,ix,0)
! CF_sum=CF_sum+2*tao0(ix)/Re
! enddo
! close(2)
! 142 format('variables = ',3(1x,A)/'zone, I =',I3,', F=POINT')
! CF_sum=CF_sum/NX
! write(*,*) NNN,CF_sum

! do iix=1,4
! write(step2,'(i2)') iix
! ix=STI(iix)
! ! open(1,file='aver1'//STEP//'.plt')
! open(1,file='aver1'//STEP//'ST'//step2//'.plt')
! write(1,131) '"y" ,','"averu" ,','"prms" ,','"Urms" ,','"Vrms" ,','"Wrms" ,','"ShearF" ,','"UV" ,',NY+1
! do iy=0,NY
! write(1,*) yy(iy),aver1(1,ix,iy),prms(ix,iy),urms(ix,iy),vrms(ix,iy),wrms(ix,iy),shearF(ix,iy),-uv(ix,iy) !
! enddo
! 131 format('variables = ',8(1x,A)/'zone, I =',I3,', F=POINT')
! close(1)
! 
! open(8,file='Uprofile'//STEP//'ST'//step2//'.plt')
! write(8,138) '"y" ,','"U+" ',NY/2+1
! do iy=0,NY/2
! write(8,*) yp(ix,iy),averu(ix,iy)
! enddo
! 138 format('variables = ',2(1x,A)/'zone, I =',I3,', F=POINT')
! close(8)
! 
! ! open(2,file='nearwall'//STEP//'.plt')
! open(2,file='nearwall'//STEP//'ST'//step2//'.plt')
! write(2,132) '"y" ,','"U+" ,','"Urms" ,','"Vrms" ,','"Wrms" ,','"ShearF" ,','"UV" ,',NY/2+1
! do iy=0,NY/2
! write(2,*) yy(iy),averu(ix,iy),urms(ix,iy),vrms(ix,iy),wrms(ix,iy),shearF(ix,iy),-uv(ix,iy)
! enddo
! 132 format('variables = ',7(1x,A)/'zone, I =',I3,', F=POINT')
! close(2)
! 
! ! open(3,file='uu'//STEP//'.plt')
! open(3,file='uu'//STEP//'ST'//step2//'.plt')
! write(3,133) '"y" ,','"PR" ,','"TD" ,','"VD" ,','"VPs" ,','"DS" ,',NY/2+1  !'"VPd" ,',
! do iy=0,NY/2
! write(3,*) yy(iy),P_11(ix,iy),T_11(ix,iy),D_11(ix,iy),PS_11(ix,iy),ipsron_11(ix,iy) !PS_11(ix,iy)+PD_11(ix,iy)
! enddo
! close(3)
! 
! ! open(4,file='vv'//STEP//'.plt')
! open(4,file='vv'//STEP//'ST'//step2//'.plt')
! write(4,133) '"y" ,','"PR" ,','"TD" ,','"VD" ,','"VPs" ,','"DS" ',NY/2+1 !,'"VPd" '
! do iy=0,NY/2
! write(4,*) yy(iy),P_22(ix,iy),T_22(ix,iy),D_22(ix,iy),PS_22(ix,iy),ipsron_22(ix,iy) !+PD_22(ix,iy)
! enddo
! ! 134 format('variables = ',6(1x,A)/'zone, I =',I3,', F=POINT')
! close(4)
! 
! ! open(5,file='ww'//STEP//'.plt')
! open(5,file='ww'//STEP//'ST'//step2//'.plt')
! write(5,133) '"y" ,','"PR" ,','"TD" ,','"VD" ,','"VPs" ,','"DS" ',NY/2+1 !,'"VPd" '
! do iy=0,NY/2
! write(5,*) yy(iy),P_33(ix,iy),T_33(ix,iy),D_33(ix,iy),PS_33(ix,iy),ipsron_33(ix,iy) !+PD_33(ix,iy)
! enddo
! ! 135 format('variables = ',5(1x,A)/'zone, I =',I3,', F=POINT')
! close(5)
! 
! ! open(6,file='uv'//STEP//'.plt')
! open(6,file='uv'//STEP//'ST'//step2//'.plt')
! write(6,133) '"y" ,','"PR" ,','"TD" ,','"VD" ,','"VPs" ,','"DS" ',NY/2+1 !,'"VPd" ,'
! do iy=0,NY/2
! write(6,*) yy(iy),P_12(ix,iy),T_12(ix,iy),D_12(ix,iy),PS_12(ix,iy),ipsron_12(ix,iy) !+PD_12(ix,iy)
! enddo
! ! 136 format('variables = ',6(1x,A)/'zone, I =',I3,', F=POINT')
! close(6)
! 
! ! open(7,file='kt'//STEP//'.plt')
! open(7,file='kt'//STEP//'ST'//step2//'.plt')
! write(7,133) '"y" ,','"PR" ,','"TD" ,','"VD" ,','"VPs" ,','"DS" ',NY/2+1 !,'"VPd" ,'
! do iy=0,NY/2
! write(7,*) yy(iy),P_k(ix,iy),T_k(ix,iy),D_k(ix,iy),PD_k(ix,iy),ipsron_k(ix,iy) !+PD_12(ix,iy)
! enddo
! close(7)
! 133 format('variables = ',6(1x,A)/'zone, I =',I3,', F=POINT')
! enddo

deallocate(xx,yy,yp,prms,averu,shearF,uu)
deallocate(urms,vrms,wrms,uv)
! deallocate(P_11,T_11,D_11,PS_11,PD_11,ipsron_11)
! deallocate(P_22,T_22,D_22,PS_22,PD_22,ipsron_22)
! deallocate(P_33,T_33,D_33,PS_33,PD_33,ipsron_33)
! deallocate(P_12,T_12,D_12,PS_12,PD_12,ipsron_12)
! deallocate(P_k ,T_k ,D_k ,PD_k ,ipsron_k)
end
