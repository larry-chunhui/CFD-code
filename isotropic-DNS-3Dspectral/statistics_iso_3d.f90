! ---&---1---------2---------3---------4---------5---------6---------7--
!计算能谱的方法5.15 周春锋
SUBROUTINE CAL_EK(U)                   !计算能量在谱空间的密度EK(I)
USE FFT
USE STATS_SPECTRUM
include 'mpif.h'
integer status(MPI_Status_Size)
integer ierr

INTEGER::KX,KY,KZ,I
COMPLEX*16 U(1:3,0:NX/2,0:NY-1,0:NZP-1)

REAL*8 TMP(0:NX/2,0:NY-1,0:NZP-1)     !计算能谱专用存储
REAL*8 EP

INTEGER COUNT            !并行

EKP=0.0D0
EK=0.0D0  !
DO KX=0,NX/2
DO KY=0,NY-1
DO KZ=0,NZP-1

   I=NINT(SQRT(real(WAVX(KX)**2+WAVY(KY)**2+WAVZ(KZ)**2)))

   IF(I.GE.0.AND.I.LE.KMAX) THEN
      TMP(KX,KY,KZ)=CDABS(U(1,KX,KY,KZ))**2+CDABS(U(2,KX,KY,KZ))**2+CDABS(U(3,KX,KY,KZ))**2
      IF(KX.EQ.0) TMP(KX,KY,KZ)=TMP(KX,KY,KZ)/2.0D0           !0平面减半
      EKP(I)=EKP(I)+TMP(KX,KY,KZ)
   ENDIF

ENDDO
ENDDO
ENDDO

! EK=EKP     !串行版本
!将所有进程的EK归纳到0进程中
COUNT=KMAX
CALL MPI_REDUCE(EKP,EK,COUNT,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

RETURN
END

!-------------------------------------------------------
SUBROUTINE CAL_SPECTRUM
USE FFT
USE STATS_SPECTRUM
USE VELOCITY
INCLUDE 'com_iso_3d.h'
INTEGER I

CALL CAL_EK(V) 

!在0进程上根据能谱计算其它统计量
IF(my_id.EQ.0) THEN

!计算X轴上的一维能谱
DO I=1,KMAX
   EKX(I)=ABS(V(1,I,0,0))**2
ENDDO

!计算总能量
EKT=0.0D0              
DO I=1,KMAX-1
   EKT=EKT+EK(I)
ENDDO
EKT=EKT+0.5D0*(EK(0)+EK(KMAX))

!计算湍动能耗散率
DISSIP=0.0D0           
DO I=1,KMAX
   DISSIP=DISSIP+I**(2.0D0)*EK(I)/Re*2.0D0
ENDDO

!计算泰勒微尺度
TaylorL=DSQRT(10.0D0*EKT/Re/DISSIP)  

!泰勒微尺度下的湍流雷诺数
Retaylor=TaylorL*DSQRT(2.0D0/3.0D0*EKT)*Re


!计算Kolmogorov尺度
KolmL=(Re**3*DISSIP)**(-1.0D0/4.0D0)


!----------------
!积分尺度和涡周转时间有问题
!----------------
!计算积分尺度,CHANGE ON 6.22
! IntegrL=(2.0D0/3.0D0*EKT)**(3.0D0/2.0D0)/DISSIP
IntegrL=0.0D0           
DO I=1,KMAX-1
   IntegrL=IntegrL+EK(I)/DFLOAT(I)           !K(I)*(K(I+1)-K(I))
ENDDO
IntegrL=IntegrL*0.75D0*PI/EKT
! 
! !计算涡周转时间
Turntime=IntegrL/(2.0D0/3.0D0*EKT)**(1.0D0/2.0D0)



!计算归一化能谱REF(张兆顺 P81)
DO I=1,KMAX
   KG(I)=DFLOAT(I)*KolmL
   EKG(I)=EK(I)*((Re)**(5.0D0/4.0D0))/((DISSIP)**(1.0D0/4.0D0))*KG(I)**(5.0D0/3.0D0)
ENDDO


ENDIF


RETURN 
END





!-------------------------------------------------------------------
!前两个波数下注入能量,维持EK(1),EK(2)不变
! REF # AMS(2015)31(1):25-31 JIN guo dong
!对并行串行通用
SUBROUTINE CAL_FORCE(V)
USE FFT
USE STATS_SPECTRUM
INCLUDE 'com_iso_3d.h'

include 'mpif.h'
integer status(MPI_Status_Size)
integer ierr

INTEGER KX,KY,KZ,I,J
REAL*8 A,B
COMPLEX*16 V(3,0:NX/2,0:NY-1,0:NZP-1)
A=0.55544D0
B=0.159843D0
EK=0.0D0

!---------------------
!计算EK(1),EK(2)，并广播给其他进程
!这一步太浪费了
!---------------------
CALL CAL_EK(V) !EK只存储在0进程中，需要BCAST给其他所有进程
CALL MPI_BCAST(EK(1),1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(EK(2),1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!----------------------
!更新K=1，2球壳内的速度
!----------------------
DO KX = 0,NX/2
DO KY = 0,NY-1
DO KZ = 0,NZP-1
    I=NINT(DSQRT(WAVX(KX)**2+WAVY(KY)**2+WAVZ(KZ)**2))
	IF(I.EQ.1) THEN
	   DO J=1,3
	      V(J,KX,KY,KZ)=DSQRT(A/EK(1))*V(J,KX,KY,KZ)
	   ENDDO
	ELSE IF(I.EQ.2) THEN
	   DO J=1,3
	      V(J,KX,KY,KZ)=DSQRT(B/EK(2))*V(J,KX,KY,KZ)
	   ENDDO
	ENDIF
ENDDO
ENDDO
ENDDO


RETURN
END
!------------------------------------------------------------------
!计算物理空间的统计量
SUBROUTINE CALP_ENERGY(VP)     !计算平均U**2,V**2,W**2
USE FFT
USE STATS_PHYSICAL
INCLUDE 'com_iso_3d.h'

include 'mpif.h'
integer status(MPI_Status_Size)
integer ierr

INTEGER I,J,L
REAL*8  Vp(3,0:NX-1,0:NYP-1,0:NZ-1)   !VP定义同谱空间相同了？

VT2P=0.0D0

DO I=0,NX-1
DO J=0,NYP-1
DO L=0,NZ-1
   VT2P(1)=VT2P(1)+VP(1,I,J,L)**2
   VT2P(2)=VT2P(2)+VP(2,I,J,L)**2
   VT2P(3)=VT2P(3)+VP(3,I,J,L)**2
ENDDO
ENDDO
ENDDO
VT2P(1)=VT2P(1)/DFLOAT(NX*NY*NZ)
VT2P(2)=VT2P(2)/DFLOAT(NX*NY*NZ)
VT2P(3)=VT2P(3)/DFLOAT(NX*NY*NZ)


!将所有进程的VT2(1),VT2(2),VT2(3)归纳到0进程中
CALL MPI_REDUCE(VT2P(1),VT2(1),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(VT2P(2),VT2(2),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(VT2P(3),VT2(3),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!CALL MPI_REDUCE(EKP,EK,COUNT,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

RETURN
END

SUBROUTINE CALP_FX(VP,VT2,FRX,GRX,HRX,GFR)     !计算两点横向相关函数F(r),G(r)
USE FFT
INCLUDE 'com_iso_3d.h'
INTEGER I,J,L,K
REAL*8  VP(1:3,0:NX-1,0:NY-1,0:NZ-1),VT2(1:3),FRX(0:NX/2),GRX(0:NX/2),HRX(0:NX/2),GFR(0:NX/2)

FRX=0.0D0
GRX=0.0D0
HRX=0.0D0
DO K=0,NX/2
   DO J=0,NYP-1
   DO L=0,NZ-1
   DO I=0,NX-1-K
      FRX(K)=FRX(K)+VP(1,I,J,L)*VP(1,I+K,J,L)
      GRX(K)=GRX(K)+VP(2,I,J,L)*VP(2,I+K,J,L)
      HRX(K)=HRX(K)+VP(3,I,J,L)*VP(3,I+K,J,L)
   ENDDO
   ENDDO
   ENDDO
   FRX(K)=FRX(K)/DFLOAT((NX-K+1)*NY*NZ)/VT2(1)
   GRX(K)=GRX(K)/DFLOAT((NX-K+1)*NY*NZ)/VT2(2)
   HRX(K)=HRX(K)/DFLOAT((NX-K+1)*NY*NZ)/VT2(3)
ENDDO


!利用f(r)计算fgr，同gr比对
DO I=1,NX/2-1
   GFR(I)=FRX(I)+I*(FRX(I+1)-FRX(I-1))/4
ENDDO
I=0
GFR(I)=FRX(I)
I=NX/2
GFR(I)=FRX(I)+DFLOAT(I)*(FRX(I)-FRX(I-1))

RETURN
END
!-----------------------------------------------------------------------
SUBROUTINE CALP_FY(VP,VT2,FRY,GRY,HRY)     !计算两点纵向相关函数FRY(r),GRY(r),HRY(r)
USE FFT
INCLUDE 'com_iso_3d.h'
INTEGER I,J,L,K
REAL*8  VP(1:3,0:NX-1,0:NY-1,0:NZ-1),VT2(1:3),FRY(0:NY/2),GRY(0:NY/2),HRY(0:NY/2),GFR(0:NX/2)

FRY=0.0D0
GRY=0.0D0
HRY=0.0D0
DO K=0,NY/2
   DO I=0,NX-1
   DO L=0,NZ-1
   DO J=0,NY-1-K
!       WRITE(*,*)'I,J,K=',I,J,K
      FRY(K)=FRY(K)+VP(1,I,J,L)*VP(1,I,J+K,L)
      GRY(K)=GRY(K)+VP(2,I,J,L)*VP(2,I,J+K,L)
      HRY(K)=HRY(K)+VP(3,I,J,L)*VP(3,I,J+K,L)
   ENDDO
   ENDDO
   ENDDO
   FRY(K)=FRY(K)/DFLOAT((NY-K+1)*NX*NZ)/VT2(1)
   GRY(K)=GRY(K)/DFLOAT((NY-K+1)*NX*NZ)/VT2(2)
   HRY(K)=HRY(K)/DFLOAT((NY-K+1)*NX*NZ)/VT2(3)
ENDDO


!利用f(r)计算fgr，同gr比对
! DO I=1,NY/2-1
!    GFR(I)=FR(I)+I*(FR(I+1)-FR(I-1))/4
! ENDDO
! I=0
! GFR(I)=FR(I)
! I=NX/2
! GFR(I)=FR(I)+I*(FR(I)-FR(I-1))

RETURN
END

!------------------------------------------------------------------------
SUBROUTINE CALP_FZ(VP,VT2,FRZ,GRZ,HRZ)     !计算两点展向Z相关函数R(i,j)
USE FFT
INCLUDE 'com_iso_3d.h'
INTEGER I,J,L,K
REAL*8  VP(1:3,0:NX-1,0:NY-1,0:NZ-1),VT2(3),FRZ(0:NZ/2),GRZ(0:NZ/2),HRZ(0:NZ/2)

FRZ=0.0D0
GRZ=0.0D0
HRZ=0.0D0
DO K=0,NZ/2
   DO I=0,NX-1
   DO J=0,NY-1
   DO L=0,NZ-1-K
      FRZ(K)=FRZ(K)+VP(1,I,J,L)*VP(1,I,J,L+K)
      GRZ(K)=GRZ(K)+VP(2,I,J,L)*VP(2,I,J,L+K)
      HRZ(K)=HRZ(K)+VP(3,I,J,L)*VP(3,I,J,L+K)
   ENDDO
   ENDDO
   ENDDO
   FRZ(K)=FRZ(K)/DFLOAT((NZ-K+1)*NX*NY)/VT2(1)
   GRZ(K)=GRZ(K)/DFLOAT((NZ-K+1)*NX*NY)/VT2(2)
   HRZ(K)=HRZ(K)/DFLOAT((NZ-K+1)*NX*NY)/VT2(3)
ENDDO

RETURN 
END


!-----------------------------------------------------------------------
SUBROUTINE CALP_Skew(VP,VT2,SKE,FLAT)                 !计算速度场偏斜因子,平坦因子
USE FFT
INCLUDE 'com_iso_3d.h'
INTEGER I,J,L,K
REAL*8  VP(1:3,0:NX-1,0:NY-1,0:NZ-1),SKE(1:4),FLAT(1:4),VT2(1:3)
SKE=0.0D0
FLAT=0.0D0

DO K=1,3
   DO I=0,NX-1
   DO J=0,NY-1
   DO L=0,NZ-1
      SKE(K)=SKE(K)+VP(K,I,J,L)**3
	  FLAT(K)=FLAT(K)+VP(K,I,J,L)**4
   ENDDO
   ENDDO
   ENDDO
   SKE(K) = SKE(K)/(VT2(k)**(1.50D0))/DFLOAT(NX*NY*NZ)
   FLAT(K)=FLAT(K)/(VT2(k)**2)/DFLOAT(NX*NY*NZ)
ENDDO
SKE(4)=(SKE(1)+SKE(2)+SKE(3))/3.0D0
FLAT(4)=(FLAT(1)+FLAT(2)+FLAT(3))/3.0D0
RETURN
END
!---------------------------------------------------------------------
SUBROUTINE CALP_DUSkew(V,SKEDU,FLATDU)                 !fft,计算速度场拉伸率的偏斜因子,平坦因子
USE FFT
INCLUDE 'com_iso_3d.h'
INTEGER I,J,L,K
COMPLEX*16 V(3,0:NX/2,0:NY-1,0:NZP-1),DV(3,0:NX/2,0:NY-1,0:NZP-1)
REAL*8 DVP(1:3,0:NX-1,0:NY-1,0:NZ-1),SKEDU(1:4),FLATDU(1:4)
REAL*8 DUT
DO I=1,3
   CALL DERIV (I,DV(I,0:NX/2,0:NY-1,0:NZP-1),V(I,0:NX/2,0:NY-1,0:NZP-1))
ENDDO
DO I=1,3
   CALL xyzfft ('B',DVP(I,0:NX-1,0:NY-1,0:NZ-1),DV(I,0:NXP/2,0:NY-1,0:NZ-1))
ENDDO
SKEDU=0.0D0
FLATDU=0.0D0
DO K=1,3
DUT=0.0D0
   DO I=0,NX-1
   DO J=0,NY-1
   DO L=0,NZ-1
      SKEDU(K)=SKEDU(K)+DVP(K,I,J,L)**3
	  FLATDU(K)=FLATDU(K)+DVP(K,I,J,L)**4
	  DUT=DUT+DVP(K,I,J,L)**2
   ENDDO
   ENDDO
   ENDDO
   DUT=DUT/DFLOAT(NX*NY*NZ)
   SKEDU(K) = SKEDU(K)/(DUT**(1.50D0))/DFLOAT(NX*NY*NZ)
   FLATDU(K)=FLATDU(K)/(DUT**2)/DFLOAT(NX*NY*NZ)
ENDDO
SKEDU(4)=(SKEDU(1)+SKEDU(2)+SKEDU(3))/3.0D0
FLATDU(4)=(FLATDU(1)+FLATDU(2)+FLATDU(3))/3.0D0
RETURN
END


! SUBROUTINE CAL_EP(EP,U)                     !计算物理空间出时空能量EP
! USE FFT
! REAL*8 EP,U(1:3,0:NX-1,0:NY-1,0:NZ-1)
! INTEGER I,J,K
! EP=0.0D0                   !COMPUTE TOTAL ENERGY 
! DO I=0,NX-1
! DO J=0,NY-1
! DO K=0,NZ-1
!    EP=EP+(U(1,I,J,K)**2+U(2,I,J,K)**2+U(3,I,J,K)**2)/2.0D0/DFLOAT(NX*NY*NZ)
! ENDDO
! ENDDO
! ENDDO
! WRITE(*,*) 'TOTAL ENERGY IN PHYSICAL SPACE'
! WRITE(*,*)EP
! 
! END

!------------------------------------------------------------------------
