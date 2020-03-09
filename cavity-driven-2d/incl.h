implicit none
integer, PARAMETER::NI=80,NJ=80
real*8,  parameter::RE=1000.0d0
integer, parameter:: Id_Escheme=1          !1:one upwind 2:hybrid 3:quick
integer, parameter:: Id_Dscheme=2          !2 order  or 4 order nearly the same!
integer, parameter::Ihalo=2,Jhalo=2

      
COMMON/VARS/U,V,P
real*8    U(-Ihalo:NI+Ihalo,-Jhalo:NJ+Jhalo),V(-Ihalo:NI+Ihalo,-Jhalo:NJ+Jhalo),P(-Ihalo:NI+Ihalo,-Jhalo:NJ+Jhalo)
COMMON/GEOM/X,Y
real*8    X(-Ihalo:NI+Ihalo,-Jhalo:NJ+Jhalo),Y(-Ihalo:NI+Ihalo,-Jhalo:NJ+Jhalo)

common/Fe/fe_u,fw_u,fn_u,fs_u,fe_v,fw_v,fn_v,fs_v
real*8 fe_u(0:1),fw_u(0:1),fn_u(0:1),fs_u(0:1)
real*8 fe_v(0:1),fw_v(0:1),fn_v(0:1),fs_v(0:1)
COMMON/COEF/APU_con,APV_con,AEU_con,AEV_con,AWU_con,AWV_con,ANU_con,ANV_con,ASU_con,ASV_con
COMMON/COEF/      AEEU_con,AEEV_con,AWWU_con,AWWV_con,ANNU_con,ANNV_con,ASSU_con,ASSV_con
real*8   APU_con(0:1),APV_con(0:1), & !0: left, 1: right for U
        AEU_con(0:1),AEV_con(0:1), &  !0: bottom, 1: top for V
        AWU_con(0:1),AWV_con(0:1), &
        ANU_con(0:1),ANV_con(0:1), &
        ASU_con(0:1),ASV_con(0:1)
real*8  AEEU_con(0:1),AEEV_con(0:1), &  !0: bottom, 1: top for V
        AWWU_con(0:1),AWWV_con(0:1), &
        ANNU_con(0:1),ANNV_con(0:1), &
        ASSU_con(0:1),ASSV_con(0:1)

COMMON/COEF/ AP_diff,AE_diff,AW_diff,AN_diff,AS_diff
COMMON/COEF/AEE_diff,AWW_diff,ANN_diff,ASS_diff
real*8 AP_diff,AE_diff,AW_diff,AN_diff,AS_diff
real*8         AEE_diff,AWW_diff,ANN_diff,ASS_diff

COMMON/COEF/APU,APV,AEU,AEV,AWU,AWV,ANU,ANV,ASU,ASV,BU,BV
COMMON/COEF/      AEEU,AEEV,AWWU,AWWV,ANNU,ANNV,ASSU,ASSV
real*8   APU(0:1),APV(0:1), & !0: left, 1: right for U
        AEU(0:1),AEV(0:1), &  !0: bottom, 1: top for V
        AWU(0:1),AWV(0:1), &
        ANU(0:1),ANV(0:1), &
        ASU(0:1),ASV(0:1), &
         BU(0:1), BV(0:1)   
real*8  AEEU(0:1),AEEV(0:1), &  !0: bottom, 1: top for V
        AWWU(0:1),AWWV(0:1), &
        ANNU(0:1),ANNV(0:1), &
        ASSU(0:1),ASSV(0:1)
COMMON/PARA_I/MAXIT
integer       MAXIT

COMMON/PARA_R/  DX,DY,URFU,URFV,URFP,SORMAX
real*8			DX,DY,URFU,URFV,URFP,SORMAX
COMMON/DEBG/RESU,RESV,RESM
real*8      RESU,RESV,RESM    

