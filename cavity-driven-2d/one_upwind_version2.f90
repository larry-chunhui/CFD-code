subroutine one_upwind
include 'incl.h'
integer k
do k=0,1
      AEU_con(k)=max(  0.0d0,-fe_u(k) )
      AWU_con(k)=max(  0.0d0, fw_u(k) )
      ANU_con(k)=max(  0.0d0,-fn_u(k) )
      ASU_con(k)=max(  0.0d0, fs_u(k) )

      AEV_con(k)=max(  0.0d0,-fe_v(k)  )
      AWV_con(k)=max(  0.0d0, fw_v(k)  )
      ANV_con(k)=max(  0.0d0,-fn_v(k)  )
      ASV_con(k)=max(  0.0d0, fs_v(k)  )

		AEEU_con(k)=0.0d0;	AEEV_con(k)=0.0d0 !0: bottom, 1: top for V
        AWWU_con(k)=0.0d0;  AWWV_con(k)=0.0d0
        ANNU_con(k)=0.0d0;  ANNV_con(k)=0.0d0
        ASSU_con(k)=0.0d0;  ASSV_con(k)=0.0d0
enddo

return
end