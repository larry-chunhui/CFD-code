! ---&---1---------2---------3---------4---------5---------6---------7--
!该版本目前支持续算，128*128*128，16进程并行
!将MODULE RS_budget改为X,Z平均
!9.22修正了Um的计算方法
!10.18添加IMPLICIT NONE
!10.23已经证明湍动能求解方案正确
!10.24搞定了脉动涡量均方根
	  MODULE FFT
	  INTEGER,SAVE :: NX,NY,NZ,NPROC,NX2,NY2,NZ2,NXP,NZP,NX2P,NZ2P,my_id,L0,NXPNY1NZP,NXPNY1NZ2P
	  INTEGER,ALLOCATABLE,SAVE :: LSUB(:)
	  REAL*8,ALLOCATABLE,SAVE :: YC(:),WAVX(:),WAVZ(:),trigx(:),trigy(:),trigz(:),trigzr(:),trigx2(:),trigy2(:),trigz2(:),WL(:),WD(:),WR(:)
	  COMPLEX*16,ALLOCATABLE,SAVE :: CWAVX(:),CWAVXRe(:),CWAVZ(:),CWAVZRe(:)
	  END MODULE

	  MODULE RS_budget 
      REAL*8,ALLOCATABLE,SAVE :: sum0(:),sum1(:,:),sum2(:,:),sum3(:,:),sum4(:,:),sump(:,:),sumdud(:,:),sumduud(:,:),sumduuud(:,:),&
                           sumddudd(:,:),sumdduudd(:,:),sumdd(:,:),sumwp(:,:),sumwp2(:,:)
      REAL*8,ALLOCATABLE,SAVE :: aver0(:),aver1(:,:),aver2(:,:),aver3(:,:),aver4(:,:),averp(:,:),averdud(:,:),averduud(:,:),averduuud(:,:),&
                           averddudd(:,:),averdduudd(:,:),averdd(:,:),averwp(:,:),averwp2(:,:)
      REAL*8,SAVE:: EKT,EKT1,EKT2,EKT3
      ! integer,save :: STI(17)
      ! integer,save :: nnx
	  END MODULE

	  MODULE plas_v
      ! REAL*8,ALLOCATABLE,SAVE :: x(:,:,:),y(:,:,:),z(:,:,:)           !网格坐标
      ! REAL*8,ALLOCATABLE,SAVE :: plas_E(:,:,:),plas_Fp(:,:,:,:)
      ! COMPLEX*16, ALLOCATABLE,SAVE :: plas_F(:,:,:,:)
	  COMPLEX*16, ALLOCATABLE :: plas_F(:,:,:,:)
	  END MODULE

      PROGRAM CHANNELDNS3D

	  USE FFT
	  use plas_v
	  use RS_budget
	  INCLUDE 'dim.h'

  include 'mpif.h'
  integer status(MPI_Status_Size)

	  INTEGER :: RSnum,NNUM0,NNUM1,NNUM2,NNUM3,NNUM4,NNUM5,NNUM6,NNUM7,NNUM8,NNUM9

	  character*2,external :: itochar



	  REAL*8, ALLOCATABLE :: Vp(:,:,:,:),VBCp(:,:,:,:)
	  REAL*8, ALLOCATABLE :: WN001(:,:),AD001(:,:,:),WNG10(:,:),ADG10(:,:,:)
	  REAL*8, ALLOCATABLE :: uvwpc(:,:,:,:,:),DVC(:,:,:,:),uG1(:)

	  COMPLEX*16, ALLOCATABLE :: V(:,:,:,:),Vm(:,:,:,:),VCW(:,:,:,:),VCWm(:,:,:,:),Wo(:,:,:,:),PRE(:,:,:),TEMP(:,:,:,:)
	  COMPLEX*16, ALLOCATABLE :: VBC(:,:,:,:),PBC(:,:,:),PBCm(:,:,:)

	  REAL*8 a(0:1),b(0:1),f
	  REAL*8 , ALLOCATABLE :: BUFFER(:)
	  INTEGER ierr

	  REAL*8 EKT_myid


! 	  REAL*8, ALLOCATABLE :: savq(:,:),RrEx(:,:,:),RrEz(:,:,:),savqtotal(:,:),RrExtotal(:,:,:),RrEztotal(:,:,:)
! 
! 	  COMPLEX*16, ALLOCATABLE :: Ex(:,:,:),Ey(:,:,:),Ez(:,:,:)

! 	  ALLOCATE(savq(49,0:NYP),RrEx(0:NXP/2-1,35,0:NYP),RrEz(0:NZ/2,35,0:NYP), &
! 			   savqtotal(49,0:NYP),RrExtotal(0:NXP/2-1,35,0:NYP),RrEztotal(0:NZ/2,35,0:NYP))
! 
! 	  ALLOCATE(Ex(0:NXP/2-1,0:NY,0:NZP-1),Ey(0:NXP/2-1,0:NY,0:NZP-1),Ez(0:NXP/2-1,0:NY,0:NZP-1))

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_Comm_World,my_id,ierr)
  call MPI_Comm_Size(MPI_Comm_World,NPROC,ierr)
!  	  my_id=0;  NPROC=1




      OPEN(10,FILE='control.dat',STATUS='UNKNOWN',ACTION='READ')
      READ(10,*)
	  READ(10,*) NX,NY,NZ
      READ(10,*)
	  READ(10,*) ZL,XL,Re,DT
      READ(10,*)
      READ(10,*) IBC
      READ(10,*)
      READ(10,*) ICONT,NSTEPS,IDUMFRQ,IPRNFRQ
      CLOSE(10)

	  IF(my_id.eq.0) then
	     WRITE(*,*)"NX,NY,NZ"
		 WRITE(*,*)NX,NY,NZ
	     WRITE(*,*)"ZL,XL,Re,DT"
		 WRITE(*,*)ZL,XL,Re,DT
		 WRITE(*,*)'IBC=',IBC
		 WRITE(*,*)"ICONT,NSTEPS,IDUMFRQ,IPRNFRQ"
		 WRITE(*,*)ICONT,NSTEPS,IDUMFRQ,IPRNFRQ
	  ENDIF
	  NX2=3*NX/2; NY2=4*NY/2;  NZ2=3*NZ/2
	  NXP=NX/NPROC; NZP=NZ/NPROC;  NX2P=NX2/NPROC;  NZ2P=NZ2/NPROC
	  NXPNY1NZP=NXP*(NY+1)*NZP;  NXPNY1NZ2P=NXP*(NY+1)*NZ2P



  ALLOCATE(BUFFER(8*NX*NY*NZ))
  call MPI_BUFFER_ATTACH(BUFFER,64*NX*NY*NZ,ierr)

	  ALLOCATE(LSUB(0:NPROC-1),YC(0:NY),WAVX(0:NXP/2-1),WAVZ(-NZ/2+1:NZ/2-1),trigx(2*NX+15),trigy(3*NY+18), &
				trigz(4*NZ+15),trigzr(4*NZ+15),trigx2(2*NX2+15),trigy2(3*NY2+18),trigz2(4*NZ2+15),WL(0:NY),WD(0:NY),WR(0:NY), &
				CWAVX(0:NXP/2-1),CWAVXRe(0:NXP/2-1),CWAVZ(-NZ/2+1:NZ/2-1),CWAVZRe(-NZ/2+1:NZ/2-1))

	  ALLOCATE(Vp(3,0:NX-1,0:NY,0:NZP-1),VBCp(3,0:NX-1,0:1,0:NZ-1))

	  ALLOCATE(WN001(0:NXP/2-1,-NZ/2+1:NZ/2-1),AD001(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
				WNG10(0:NXP/2-1,-NZ/2+1:NZ/2-1),ADG10(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
				uvwpc(4,2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),DVC(2,0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1),uG1(0:NY))

	  ALLOCATE(V(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),Vm(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
				VCW(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),VCWm(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
				Wo(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),PRE(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))

	  ALLOCATE(VBC(3,0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1),PBC(0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1),PBCm(0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1))

allocate(sum0(2),sum1(3,0:NY)  ,sum2(6,0:NY), sum3(3,0:NY),sum4(3,0:NY), sump(10,0:NY), sumdud(9,0:NY),&
         sumduud(16,0:NY), sumduuud(12,0:NY),  sumddudd(9,0:NY),&
		 sumdduudd(12,0:NY),  sumdd(12,0:NY),  sumwp(3,0:NY),sumwp2(3,0:NY))
allocate(aver0(2),aver1(3,0:NY),aver2(6,0:NY),aver3(3,0:NY),aver4(3,0:NY),averp(10,0:NY),averdud(9,0:NY),&
         averduud(16,0:NY),averduuud(12,0:NY),averddudd(9,0:NY),&
		 averdduudd(12,0:NY),averdd(12,0:NY), averwp(3,0:NY),averwp2(3,0:NY))
ALLOCATE(plas_F(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! NNX=17
! allocate(sum1(3,NNX,0:NY)  ,sum2(6,NNX,0:NY), sump(10,NNX,0:NY), sumdud(9,NNX,0:NY),&
!          sumduud(16,NNX,0:NY), sumduuud(12,NNX,0:NY),  sumddudd(9,NNX,0:NY),&
! 		 sumdduudd(12,NNX,0:NY),  sumdd(12,NNX,0:NY))
! allocate(aver1(3,NNX,0:NY),aver2(6,NNX,0:NY),averp(10,NNX,0:NY),averdud(9,NNX,0:NY),&
!          averduud(16,NNX,0:NY),averduuud(12,NNX,0:NY),averddudd(9,NNX,0:NY),&
! 		 averdduudd(12,NNX,0:NY),averdd(12,NNX,0:NY))

	  CALL CPU_TIME(time1)

	  CALL SETUP

	  IF(ICONT==0) THEN
		 NT=0; 	TIME=0.0D0
		 CALL VBOUND (VBCp,VBC,TIME)
	     CALL INITIAL(Vp,VBCp)
	  ELSE
		 OPEN(50,FILE='field='//itochar(my_id)//'.dat',STATUS='UNKNOWN')
		 READ(50,*) NT,TIME
		 READ(50,*) Vp(1:3,0:NX-1,0:NY,0:NZP-1)      !CHANGE BY LIU CHUNHUI ON 9.16
		 CLOSE(50)
		 IF(IBC==0) CALL VBOUND (VBCp,VBC,TIME)
	  ENDIF

      DO I=1,3
		 CALL xyzfft ('F',Vp(I,0:NX-1,0:NY,0:NZP-1),V(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
      ENDDO

	  CALL SOLVE_i (WN001,AD001,0.0D0,0.0D0,1.0D0)

! 	  NSAMPL = 0
! 	  savqtotal = 0.0D0
! 	  RrExtotal = 0.0D0
! 	  RrEztotal = 0.0D0
EKT1=0.0D0
EKT2=0.0D0
EKT3=0.0D0
sum0=0.0d0
sum1=0.0D0
sum2=0.0D0
sum3=0.0d0
sum4=0.0d0
sump=0.0D0
sumdud=0.0D0
sumduud=0.0D0
sumduuud=0.0D0
sumddudd=0.0D0
sumdduudd=0.0D0
sumdd=0.0D0
sumwp=0.0d0
sumwp2=0.0d0
RSnum=0
NNUM0=2
NNUM1=3*(NY+1)
NNUM2=6*(NY+1)
NNUM3=10*(NY+1)
NNUM4=9*(NY+1)
NNUM5=16*(NY+1)
NNUM6=12*(NY+1)
NNUM7=9*(NY+1)
NNUM8=12*(NY+1)
NNUM9=12*(NY+1)
!*************************************
!       call plas_actuator
!*************************************

! call lamda(V)


!**************************
!by huangjunji  delete by liuchunhui
!**************************
      IF(my_id.eq.0) then
 		OPEN(61,FILE='EKT.plt',STATUS='UNKNOWN')
		WRITE(61,*)"VARIABLES=TIME,EKT0,EKT1,EKT2,EKT3"
	  ENDIF

      DO 290 it=1,NSTEPS
	     IF(it<3) THEN
			a(0) = DFLOAT(it);  a(1) = DFLOAT(1-it)/2.0D0
			b(0) = DFLOAT(it);  b(1) = DFLOAT(1-it)
			a = a/DT       !已经除以了一个DT
			G = (a(0)+a(1))*Re
			CALL Q_c (uG1,QG1,G)
			CALL DIV_c (uvwpc,DVC,G,WNG10,ADG10)
		 ENDIF

         TIME = TIME + DT
		 IF(my_id .EQ. 0) WRITE(*,*) 'it = ', it,'    TIME = ', TIME

		 IF(IBC==1) CALL VBOUND (VBCp,VBC,TIME)

! 		 DO I=1,3
! 			DO J=1,3
! 			   CALL DERIV (J,gradV(I,J,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),V(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
! 			ENDDO
! 		 ENDDO
		 DO I=1,3     !求涡量
			K1=E_ijk(I,1);  K2=E_ijk(I,2)
			CALL DERIV (K1,TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),V(K2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
			CALL DERIV (K2,TEMP(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),V(K1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
			Wo(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) = TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) &
												- TEMP(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)
! 			Wo(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) = gradV(K1,K2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) &
! 												- gradV(K2,K1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)
		 ENDDO

		 CALL CROSSTIME(VCW,V,Wo)    !非线性项，ucw=u*w

		 TEMP = (a(0)*V+a(1)*Vm) + (b(0)*VCW+b(1)*VCWm)
! 		 TEMP = (a(0)*V+a(1)*Vm) + (b(0)*VCW+b(1)*VCWm) + plas_F

		 Vm=V;  VCWm=VCW;  V=TEMP

	     PBC = DCMPLX(0.0D0, 0.0D0)
	     DO iy=0,NY      !压力边界条件，PBC(0)=a+,y=1(iy=0);PBC(1)=a-,y=-1(iy=NY)。
			DO iz=-NZ/2+1,NZ/2-1
			   DO ix=0,NXP/2-1
				  TEMP(1,ix,0,iz)=CWAVXRe(ix)*Wo(3,ix,iy,iz)-CWAVZRe(iz)*Wo(1,ix,iy,iz)+VCW(2,ix,iy,iz)
			   ENDDO
			ENDDO
			PBC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)=PBC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)+TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)
			IF(MOD(iy,2)==0) THEN
			   PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)=PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)-TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)
			ELSE
			   PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)=PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)+TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)
			ENDIF
	     ENDDO

		 DO I=1,3
			CALL DERIV (I,TEMP(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),V(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
		 ENDDO
	     TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) = TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) &
												+TEMP(2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) &
												+TEMP(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)

	     PBCm = b(0)*PBC+b(1)*PBCm
	     CALL SOLVE (PRE,TEMP(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),0.0D0,0.0D0,1.0D0,PBCm,WN001,AD001)
	     PBCm = PBC

		 DO I=1,3
			CALL DERIV (I,TEMP(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),PRE)
		 ENDDO
         TEMP = Re*(TEMP-V)

		 DO I=1,3        !temp为Helmholtz方程的右端项，WNG10
			CALL SOLVE(V(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),TEMP(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
									G,1.0D0,0.0D0,VBC(I,0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1),WNG10,ADG10)
	     ENDDO
	     PBC = DCMPLX(0.0D0, 0.0D0)
	     DO iy=0,NY
	        TEMP(1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)=DFLOAT(iy**2)*V(2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
	        PBC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)=PBC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)+TEMP(1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
			IF(MOD(iy,2)==0) THEN
			   PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)=PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)-TEMP(1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
			ELSE
	           PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)=PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)+TEMP(1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
			ENDIF
	     ENDDO

	     TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1) = PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)*DVC(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1) &
											- PBC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*DVC(2,0:NXP/2-1,1,-NZ/2+1:NZ/2-1)
	     TEMP(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1) = PBC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*DVC(1,0:NXP/2-1,1,-NZ/2+1:NZ/2-1) &
											- PBC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)*DVC(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)

	     DO iy = 0,NY
	        V(1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = V(1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) + DCMPLX(0.0D0,1.0D0) &
									*(TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(1,1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) &
									+ TEMP(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(1,2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1))
	        V(2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = V(2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) &
									+ TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(2,1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) &
									+ TEMP(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(2,2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
	        V(3,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = V(3,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) + DCMPLX(0.0D0,1.0D0) &
									*(TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(3,1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) &
									+ TEMP(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(3,2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1))
	        PRE(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = PRE(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) &
									+ TEMP(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(4,1,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) &
									+ TEMP(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*uvwpc(4,2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
	     ENDDO

!###################################################
!uG1是f=1时对应的流速，QG1是对应的流量，都不随时间变化，可一次性求出。
!流量Q,UG1,F是物理空间的量，只是计算流量时用到谱空间V(1,0,0:NY,0)来计算
!看起来在谱空间进行的流量修正
!###################################################
	     IF(my_id .EQ. 0) THEN   
	        Q=0.0D0
	        DO iy=0,NY,2
	           Q = Q - DREAL(V(1,0,iy,0))/DFLOAT(iy**2-1)
	        ENDDO
	        f = (1.0D0-Q)/QG1
	        V(1,0,0:NY,0) = V(1,0,0:NY,0) + f*uG1
	     ENDIF



!###################################################
! if (it .ge. 2000 .and. mod(it,10) .eq. 0) then
! if (it .GE.2000.AND.MOD(it,10).eq.0 ) then
if (it .EQ.1 ) then
RSnum=RSnum+1
call RS_sum(V,PRE)
endif

! if (MOD(it,IDUMFRQ) .EQ. 0) then
if(it.eq.1) then
aver0=0.0d0
aver1=0.0d0
aver2=0.0d0
aver3=0.0d0
aver4=0.0d0
averp=0.0d0
averdud=0.0d0
averduud=0.0d0
averduuud=0.0d0
averddudd=0.0d0
averdduudd=0.0d0
averdd=0.0d0
averwp=0.0d0
averwp2=0.0d0


CALL MPI_REDUCE(sum0,aver0,NNUM0,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sum1,aver1,NNUM1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sum2,aver2,NNUM2,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sum3,aver3,NNUM1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sum4,aver4,NNUM1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sump,averp,NNUM3,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumdud,averdud,NNUM4,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumduud,averduud,NNUM5,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumduuud,averduuud,NNUM6,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumddudd,averddudd,NNUM7,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumdduudd,averdduudd,NNUM8,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumdd,averdd,NNUM9,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

CALL MPI_REDUCE(sumwp,averwp,NNUM1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(sumwp2,averwp2,NNUM1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
if (my_id .EQ. 0) THEN
aver0(1)=aver0(1)/DFLOAT(RSnum*NZ*NX)          !sum0(1)=sum0(1)+up(1,ix,NY/2,iz)    !Uc
aver0(2)=aver0(2)/DFLOAT(RSnum*NZ*NX*2)        !sum0(2)=sum0(2)+up(1,ix,iy,iz)      !Um
aver1=aver1/DFLOAT(RSnum*NZ*NX)
aver2=aver2/DFLOAT(RSnum*NZ*NX)
aver3=aver3/DFLOAT(RSnum*NZ*NX)            !<u**3>
aver4=aver4/DFLOAT(RSnum*NZ*NX)            !<u**4>
averwp=averwp/DFLOAT(RSnum*NZ*NX)          !<w>
averwp2=averwp2/DFLOAT(RSnum*NZ*NX)        !<w*w>

averp=averp/DFLOAT(RSnum*NZ*NX)
averdud=averdud/DFLOAT(RSnum*NZ*NX)
averduud=averduud/DFLOAT(RSnum*NZ*NX)
averduuud=averduuud/DFLOAT(RSnum*NZ*NX)
averddudd=averddudd/DFLOAT(RSnum*NZ*NX)
averdduudd=averdduudd/DFLOAT(RSnum*NZ*NX)
averdd=averdd/DFLOAT(RSnum*NZ*NX)
call RS_statis(RSnum)
WRITE(*,*)'Q=',Q,'f=',f
endif

endif

CALL cal_k(V)

IF(my_id.eq.0) THEN
	WRITE(61,*)TIME,EKT1+EKT2+EKT3,EKT1,EKT2,EKT3
	WRITE(*,*)EKT1,EKT2,EKT3
ENDIF
!CHECK EKT BY HUANGJUN JI
!已经证明湍动能求解方案正确

			EKT_myid=0.0D0
			EKT=0.0D0
			DO iy=0,NY
			   DO ix=0,NX-1
				  DO iz=0,NZP-1
					 DO I=1,3
						EKT_myid=EKT_myid+Vp(I,ix,iy,iz)**2
					 ENDDO
				  ENDDO
			   ENDDO
			ENDDO

			CALL MPI_REDUCE(EKT_myid,EKT,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
			IF(my_id.eq.0) then
				WRITE(*,*)"CHECK EKT BY HUANGJUN JI"
				WRITE(*,*)"EKT=",EKT/dfloat(NX*(NY+1)*NZ)
			endif
!stop !add by liu chunhui
!###################################################

!          IF(MOD(it,IDUMFRQ) .EQ. 0) THEN
           IF(it.eq.IDUMFRQ) THEN
            call lamda(V)

			DO I=1,3
			   CALL xyzfft ('B',Vp(I,0:NX-1,0:NY,0:NZP-1),V(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
			ENDDO

			OPEN(60,FILE='field='//itochar(my_id)//'.dat',STATUS='UNKNOWN')
			WRITE(60,*) NT+it,TIME
      		WRITE(60,*) Vp(1:3,0:NX-1,0:NY,0:NZP-1)
			CLOSE(60)
!**************************************************************
!disturb the flow field to reach fully developed turbulence
!**************************************************************
 			IF(NT+it<5001) THEN
 			   CALL DISTURB (Vp)
 			   DO I=1,3
 				  CALL xyzfft ('F',Vp(I,0:NX-1,0:NY,0:NZP-1),V(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
 			   ENDDO
 			ENDIF

		 ENDIF



 290  CONTINUE

	  CLOSE(61)

!       IF(NSAMPL.GT.0) THEN
!  
! 		 savq = savqtotal/DFLOAT(NSAMPL)
! 		 RrEx = RrExtotal/DFLOAT(NSAMPL)
! 		 RrEz = RrEztotal/DFLOAT(NSAMPL)
! 
! 		 DO iy=0,NY
! 			OPEN(60,FILE='avq_iy='//itochar(my_id*NPROC+iy)//'.dat',STATUS='UNKNOWN')
! 			WRITE(60,*) savq(1:49,iy),RrEx(0:NXP/2-1,1:35,iy),RrEz(0:NXP/2-1,1:35,iy)
! 			CLOSE(60)
! 		 ENDDO
! 
! 	  ENDIF
			
	  CALL CPU_TIME(time2)
      WRITE(*,*) 'Time = ', time2-time1

      call MPI_finalize(ierr)
	  
	  DEALLOCATE(LSUB,YC,WAVX,WAVZ,trigx,trigy,trigz,trigx2,trigy2,trigz2,WL,WD,WR,CWAVX,CWAVXRe,CWAVZ,CWAVZRe)

	  DEALLOCATE(Vp,VBCp,WN001,AD001,WNG10,ADG10,uvwpc,DVC,uG1)

	  DEALLOCATE(V,Vm,VCW,VCWm,Wo,PRE,TEMP)

	  DEALLOCATE(VBC,PBC,PBCm)

	  DEALLOCATE(sum1,sum2,sump,sumdud,sumduud,sumduuud,sumddudd,sumdduudd,sumdd)
      DEALLOCATE(aver1,aver2,averp,averdud,averduud,averduuud,averddudd,averdduudd,averdd)

! 	  DEALLOCATE(uvwpc,DUC,uG0,savq,RrEx,RrEz,savqtotal,RrExtotal,RrEztotal,BUFFER)
! 
! 	  DEALLOCATE(Ex,Ey,Ez)

      END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  CHARACTER*8 FUNCTION Ftochar(Q)

	  IMPLICIT REAL*8 (A-H, O-Z)

	  WRITE(Ftochar,'(F8.2)') Q

	  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  CHARACTER*2 FUNCTION itochar(iy)

	  IMPLICIT REAL*8 (A-H, O-Z)

	  WRITE(itochar,'(I2)') iy

	  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE CROSSTIME(ACB,A,B)

	  USE FFT
	  INCLUDE 'dim.h'

	  COMPLEX*16, DIMENSION(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) :: ACB,A,B

	  REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ACB2(:,:,:,:)

	  ALLOCATE(A2(3,0:NX2-1,0:NY2,0:NZ2P-1),B2(3,0:NX2-1,0:NY2,0:NZ2P-1),ACB2(3,0:NX2-1,0:NY2,0:NZ2P-1))


	  DO I=1,3
		 CALL xyzfft2 ('B',A2(I,0:NX2-1,0:NY2,0:NZ2P-1),A(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
		 CALL xyzfft2 ('B',B2(I,0:NX2-1,0:NY2,0:NZ2P-1),B(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
	  ENDDO

	  DO I=1,3
		 K1=E_ijk(I,1);  K2=E_ijk(I,2)
		 ACB2(I,0:NX2-1,0:NY2,0:NZ2P-1) = A2(K1,0:NX2-1,0:NY2,0:NZ2P-1)*B2(K2,0:NX2-1,0:NY2,0:NZ2P-1) - &
										  A2(K2,0:NX2-1,0:NY2,0:NZ2P-1)*B2(K1,0:NX2-1,0:NY2,0:NZ2P-1)
	  ENDDO

	  DO I=1,3
		 CALL xyzfft2 ('F',ACB2(I,0:NX2-1,0:NY2,0:NZ2P-1),ACB(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
	  ENDDO

	  DEALLOCATE(A2,B2,ACB2)

	  RETURN
	  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE DOTTIME(ADB,A,B)

	  USE FFT
	  INCLUDE 'dim.h'

	  COMPLEX*16 ADB(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)
	  COMPLEX*16, DIMENSION(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) :: A,B

	  REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ADB2(:,:,:)

	  ALLOCATE(A2(3,0:NX2-1,0:NY2,0:NZ2P-1),B2(3,0:NX2-1,0:NY2,0:NZ2P-1),ADB2(0:NX2-1,0:NY2,0:NZ2P-1))

	  DO I=1,3
		 CALL xyzfft2 ('B',A2(I,0:NX2-1,0:NY2,0:NZ2P-1),A(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
		 CALL xyzfft2 ('B',B2(I,0:NX2-1,0:NY2,0:NZ2P-1),B(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
	  ENDDO

	  ADB2 = A2(1,0:NX2-1,0:NY2,0:NZ2P-1)*B2(1,0:NX2-1,0:NY2,0:NZ2P-1) &
			+A2(2,0:NX2-1,0:NY2,0:NZ2P-1)*B2(2,0:NX2-1,0:NY2,0:NZ2P-1) &
			+A2(3,0:NX2-1,0:NY2,0:NZ2P-1)*B2(3,0:NX2-1,0:NY2,0:NZ2P-1)


	  CALL xyzfft2 ('F',ADB2(0:NX2-1,0:NY2,0:NZ2P-1),ADB(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))

	  DEALLOCATE(A2,B2,ADB2)

	  RETURN
	  END

