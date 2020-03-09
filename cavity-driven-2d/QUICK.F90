subroutine quick_upwind
include 'incl.h'
real,external ::alfa
integer k
!test function alfa
!write(*,*)"alfa(1.00)=1??",alfa(100.0d0)
!write(*,*)"alfa(-1.00)=1??",alfa(-9.0d0)
!stop
do k=0,1
      AWU_con(k) =6.0d0/8.0d0*alfa(fw_u(k))*fw_u(k)+1.0d0/8.0d0*alfa(fe_u(k))*fe_u(k) &
					  +3.0d0/8.0d0*(1-alfa(fw_u(k)))*fw_u(k)
      AWWU_con(k)=-1.0d0/8.0d0*alfa(fw_u(k))*fw_u(k)

      AEU_con(k) =-3.0d0/8.0d0*alfa(fe_u(k))*fe_u(k)-6.0d0/8.0d0*( 1-alfa(fe_u(k)) )*fe_u(k) &
				     -1.0d0/8.0d0*( 1-alfa(fw_u(k)) )*fw_u(k)
      AEEU_con(k)=1.0d0/8.0d0*(1-alfa(fe_u(k)))*fe_u(k)

	  ASU_con(k) =6.0d0/8.0d0*alfa(fs_u(k))*fs_u(k)+1.0d0/8.0d0*alfa(fn_u(k))*fn_u(k) &
					  +3.0d0/8.0d0*(1-alfa(fs_u(k)))*fs_u(k)
	  ASSU_con(k)=-1.0d0/8.0d0*alfa(fs_u(k))*fs_u(k)
!
	  ANU_con(k)=-3.0D0/8.0D0*alfa(fn_u(k))*fn_u(k)-6.0d0/8.0d0*( 1.0d0-alfa(fn_u(k)) )*fn_u(k) &
					   -1.0d0/8.0d0*( 1.0d0-alfa(fs_u(k)) )*fs_u(k)
	  ANNU_con(k)=1.0D0/8.0D0*( 1-alfa(fn_u(k)) )*fn_u(k)


!--------------------------------------------------------------------------------------------------------
      AWV_con(k) =6.0d0/8.0d0*alfa(fw_v(k))*fw_v(k)+1.0d0/8.0d0*alfa(fe_v(k))*fe_v(k) &
					  +3.0d0/8.0d0*(1-alfa(fw_v(k)))*fw_v(k)
      AWWV_con(k)=-1.0d0/8.0d0*alfa(fw_v(k))*fw_v(k)

      AEV_con(k) =-3.0d0/8.0d0*alfa(fe_v(k))*fe_v(k)-6.0d0/8.0d0*( 1-alfa(fe_v(k)) )*fe_v(k) &
				     -1.0d0/8.0d0*( 1-alfa(fw_v(k)) )*fw_v(k)
      AEEV_con(k)=1.0d0/8.0d0*(1-alfa(fe_v(k)))*fe_v(k)

	  ASV_con(k) =6.0d0/8.0d0*alfa(fs_v(k))*fs_v(k)+1.0d0/8.0d0*alfa(fn_v(k))*fn_v(k) &
					  +3.0d0/8.0d0*(1-alfa(fs_v(k)))*fs_v(k)
	  ASSV_con(k)=-1.0d0/8.0d0*alfa(fs_v(k))*fs_v(k)

	  ANV_con(k) =-3.0D0/8.0D0*alfa(fn_v(k))*fn_v(k)-6.0d0/8.0d0*( 1.0d0-alfa(fn_v(k)) )*fn_v(k) &
					   -1.0d0/8.0d0*( 1.0d0-alfa(fs_v(k)) )*fs_v(k)
	  ANNV_con(k)=1.0D0/8.0D0*( 1-alfa(fn_v(k)) )*fn_v(k)


enddo


return
end

function alfa(a)
implicit none
real*8 :: a,alfa

if(a.gt.0.0) alfa=1.0d0
if(a.le.0.0) alfa=0.0d0
return
end