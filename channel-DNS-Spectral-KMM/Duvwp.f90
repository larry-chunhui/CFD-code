! call dxyz(up,pp,uu,dpd1,duud,duuud,dduudd)
subroutine dxyz(up,pp,uu,dpd,duud,duuud,dduudd)
use FFT
INCLUDE 'dim.h'
real*8,intent(in) :: up(3,0:NX-1,0:NY,0:NZP-1),pp(0:NX-1,0:NY,0:NZP-1),uu(6,0:NX-1,0:NY,0:NZP-1)
real*8,intent(out) :: dpd(5,0:NX-1,0:NY,0:NZP-1),duud(16,0:NX-1,0:NY,0:NZP-1),&
                      duuud(12,0:NX-1,0:NY,0:NZP-1),dduudd(12,0:NX-1,0:NY,0:NZP-1)
real*8,allocatable :: uuu(:,:,:,:),pu(:,:,:,:)
complex*16,allocatable :: temp1(:,:,:,:),temp2(:,:,:,:),temp3(:,:,:,:),temp11(:,:,:,:),temp22(:,:,:,:),temp33(:,:,:,:),temp222(:,:,:,:)
allocate(uuu(10,0:NX-1,0:NY,0:NZP-1),pu(3,0:NX-1,0:NY,0:NZP-1))
allocate(temp1(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp3(10,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),&
         temp11(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp22(16,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp33(12,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),&
		 temp222(12,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))

DO iy=0,NY
DO ix=0,NX-1
DO iz=0,NZP-1
pu(1,ix,iy,iz)=pp(ix,iy,iz)*up(1,ix,iy,iz)    !pu
pu(2,ix,iy,iz)=pp(ix,iy,iz)*up(2,ix,iy,iz)    !pv
pu(3,ix,iy,iz)=pp(ix,iy,iz)*up(3,ix,iy,iz)    !pw
uuu(1,ix,iy,iz) =uu(1,ix,iy,iz)*up(1,ix,iy,iz)    !uuu
uuu(2,ix,iy,iz) =uu(1,ix,iy,iz)*up(2,ix,iy,iz)    !uuv
uuu(3,ix,iy,iz) =uu(1,ix,iy,iz)*up(3,ix,iy,iz)    !uuw
uuu(4,ix,iy,iz) =uu(2,ix,iy,iz)*up(1,ix,iy,iz)    !vvu
uuu(5,ix,iy,iz) =uu(2,ix,iy,iz)*up(2,ix,iy,iz)    !vvv
uuu(6,ix,iy,iz) =uu(2,ix,iy,iz)*up(3,ix,iy,iz)    !vvw
uuu(7,ix,iy,iz) =uu(3,ix,iy,iz)*up(1,ix,iy,iz)    !wwu
uuu(8,ix,iy,iz) =uu(3,ix,iy,iz)*up(2,ix,iy,iz)    !wwv
uuu(9,ix,iy,iz) =uu(3,ix,iy,iz)*up(3,ix,iy,iz)    !www
uuu(10,ix,iy,iz)=uu(4,ix,iy,iz)*up(3,ix,iy,iz)    !uvw
enddo
enddo
enddo
do i =1,3
CALL xyzfft ('F',pu(i,0:NX-1,0:NY,0:NZP-1),temp1(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !pu/v/w
enddo
do i =1,6
CALL xyzfft ('F',uu(i,0:NX-1,0:NY,0:NZP-1),temp2(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !uu/v/w
enddo
do i =1,10
CALL xyzfft ('F',uuu(i,0:NX-1,0:NY,0:NZP-1),temp3(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)) !uuu/v/w
enddo

CALL DERIV (1,TEMP11(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp1(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !dpud
CALL DERIV (2,TEMP11(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp1(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP11(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp1(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP11(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp1(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP11(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp1(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i=1,5
call xyzfft ('B',dpd(i,0:NX-1,0:NY,0:NZP-1),temp11(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
enddo

CALL DERIV (1,TEMP22(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !duud
CALL DERIV (2,TEMP22(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP22(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP22(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP22(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP22(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP22(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP22(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP22(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp2(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP22(10,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP22(11,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP22(12,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP22(13,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP22(14,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP22(15,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP22(16,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp2(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i=1,16
call xyzfft ('B',duud(i,0:NX-1,0:NY,0:NZP-1),temp22(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
enddo
CALL DERIV (1,TEMP222(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !dduudd
CALL DERIV (2,TEMP222(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP222(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP222(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP222(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP222(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP222(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP222(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP222(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp22(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP222(10,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp22(10,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP222(11,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp22(11,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP222(12,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp22(12,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i=1,12
call xyzfft ('B',dduudd(i,0:NX-1,0:NY,0:NZP-1),temp222(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
enddo

CALL DERIV (1,TEMP33(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))  !duuud
CALL DERIV (2,TEMP33(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP33(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP33(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP33(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(5,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP33(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(6,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP33(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(7,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP33(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(8,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP33(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) ,temp3(9,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (1,TEMP33(10,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp3(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (2,TEMP33(11,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp3(4,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
CALL DERIV (3,TEMP33(12,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),temp3(10,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
do i=1,12
call xyzfft ('B',duuud(i,0:NX-1,0:NY,0:NZP-1),temp33(i,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
enddo

deallocate(uuu,pu,temp1,temp2,temp3,temp11,temp22,temp33,temp222)
end