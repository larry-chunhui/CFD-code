subroutine quick_upwind(I,J,b1,b2,b3,b4,b5)
include 'incl.h'
real*8 b1,b2,b3,b4,b5
real,external ::alfa
integer I,J,k
!test function alfa
!write(*,*)"alfa(1.00)=1??",alfa(100.0d0)
!write(*,*)"alfa(-1.00)=1??",alfa(-9.0d0)
!stop
do k=0,1
      AWU(k) =1.0d0/RE+6.0d0/8.0d0*alfa(fw_u(k))*fw_u(k)+1.0d0/8.0d0*alfa(fe_u(k))*fe_u(k) &
					  +3.0d0/8.0d0*(1-alfa(fw_u(k)))*fw_u(k)
      AWWU(k)=-1.0d0/8.0d0*alfa(fw_u(k))*fw_u(k)

      AEU(k) =1.0d0/RE-3.0d0/8.0d0*alfa(fe_u(k))*fe_u(k)-6.0d0/8.0d0*( 1-alfa(fe_u(k)) )*fe_u(k) &
				     -1.0d0/8.0d0*( 1-alfa(fw_u(k)) )*fw_u(k)
      AEEU(k)=1.0d0/8.0d0*(1-alfa(fe_u(k)))*fe_u(k)

	  ASU(k) =1.0d0/RE+6.0d0/8.0d0*alfa(fs_u(k))*fs_u(k)+1.0d0/8.0d0*alfa(fn_u(k))*fn_u(k) &
					  +3.0d0/8.0d0*(1-alfa(fs_u(k)))*fs_u(k)
	  ASSU(k)=-1.0d0/8.0d0*alfa(fs_u(k))*fs_u(k)
!
	  ANU(k) =1.0D0/RE-3.0D0/8.0D0*alfa(fn_u(k))*fn_u(k)-6.0d0/8.0d0*( 1.0d0-alfa(fn_u(k)) )*fn_u(k) &
					   -1.0d0/8.0d0*( 1.0d0-alfa(fs_u(k)) )*fs_u(k)
	  ANNU(k)=1.0D0/8.0D0*( 1-alfa(fn_u(k)) )*fn_u(k)

      APU(k)=AEU(k)+AWU(k)+ANU(k)+ASU(k)+AEEU(k)+AWWU(k)+ANNU(k)+ASSU(k)
!--------------------------------------------------------------------------------------------------------
      AWV(k) =1.0d0/RE+6.0d0/8.0d0*alfa(fw_v(k))*fw_v(k)+1.0d0/8.0d0*alfa(fe_v(k))*fe_v(k) &
					  +3.0d0/8.0d0*(1-alfa(fw_v(k)))*fw_v(k)
      AWWV(k)=-1.0d0/8.0d0*alfa(fw_v(k))*fw_v(k)

      AEV(k) =1.0d0/RE-3.0d0/8.0d0*alfa(fe_v(k))*fe_v(k)-6.0d0/8.0d0*( 1-alfa(fe_v(k)) )*fe_v(k) &
				     -1.0d0/8.0d0*( 1-alfa(fw_v(k)) )*fw_v(k)
      AEEV(k)=1.0d0/8.0d0*(1-alfa(fe_v(k)))*fe_v(k)

	  ASV(k) =1.0d0/RE+6.0d0/8.0d0*alfa(fs_v(k))*fs_v(k)+1.0d0/8.0d0*alfa(fn_v(k))*fn_v(k) &
					  +3.0d0/8.0d0*(1-alfa(fs_v(k)))*fs_v(k)
	  ASSV(k)=-1.0d0/8.0d0*alfa(fs_v(k))*fs_v(k)

	  ANV(k) =1.0D0/RE-3.0D0/8.0D0*alfa(fn_v(k))*fn_v(k)-6.0d0/8.0d0*( 1.0d0-alfa(fn_v(k)) )*fn_v(k) &
					   -1.0d0/8.0d0*( 1.0d0-alfa(fs_v(k)) )*fs_v(k)
	  ANNV(k)=1.0D0/8.0D0*( 1-alfa(fn_v(k)) )*fn_v(k)

      APV(k)=AEV(k)+AWV(k)+ANV(k)+ASV(k)+AEEV(k)+AWWV(k)+ANNV(k)+ASSV(k)
enddo


       BU(1)=( P(I,J)-P(I+1,J) )*DY    
!c for u(i-1,j)

       BU(0)=( P(I-1,J)-P(I,J) )*DY  
!c for v(i,j  )

       BV(1)=( P(I,J)-P(I,J+1) )*DX 	
!c for v(i,j-1)

       BV(0)=( P(I,J-1)-P(I,J) )*DX 
!U(i-1,j)
      b1=-APU(0)*U(I-1,J)+AEU( 0)*U(I,  J)+AWU( 0)*U(I-2,J)+ANU( 0)*U(I-1,J+1)+ASU( 0)*U(I-1,J-1)+BU(0) &
						 +AEEU(0)*U(I+1,J)+AWWU(0)*U(I-3,J)+ANNU(0)*U(I-1,J+2)+ASSU(0)*U(I-1,J-2)
!U(i,j)   
      b2=-APU(1)*U(I,  J)+AEU( 1)*U(I+1,J) +AWU(1)*U(I-1,J) +ANU(1)*U(I,J+1)+ASU( 1)*U(I,J-1)+BU(1) &
						 +AEEU(1)*U(I+2,J)+AWWU(1)*U(I-2,J)+ANNU(1)*U(I,J+2)+ASSU(1)*U(I,J-2)
!V(i,j-1)
      b3=-APV(0)*V(I,J-1)+AEV( 0)*V(I+1,J-1)+AWV( 0)*V(I-1,J-1)+ANV( 0)*V(I,J  )+ASV( 0)*V(I,J-2)+BV(0) &
						 +AEEV(0)*V(I+2,J-1)+AWWV(0)*V(I-2,J-1)+ANNV(0)*V(I,J+1)+ASSV(0)*V(I,J-3)
!V(i,j)
      b4=-APV(1)*V(I,  J)+AEV( 1)*V(I+1,J)+AWV( 1)*V(I-1,J)+ANV( 1)*V(I,  J+1)+ASV( 1)*V(I,  J-1)+BV(1) &
						 +AEEV(1)*V(I+2,J)+AWWV(1)*V(I-2,J)+ANNV(1)*V(I,  J+2)+ASSV(1)*V(I,  J-2)
      b5=( U(I-1,J)-U(I,J) )*DY + ( V(I,J-1)-V(I,J) )*DX

return
end

function alfa(a)
implicit none
real*8 :: a,alfa

if(a.gt.0.0) alfa=1.0d0
if(a.le.0.0) alfa=0.0d0
return
end