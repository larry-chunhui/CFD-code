!���а汾
!5.19 BY ������ 
!	  CALL CAL_EKT(EK,EKT,KMAX)                     !�����ʼ�Ķ���EKT0 
!     CALL CAL_dissipation(EK,KMAX,DISSIP)          !�����Ķ��ܺ�ɢ��
!     CALL CAL_TaylorL(EKT,DISSIP,TaylorL)          !����̩��΢�߶�
!     CALL CAL_Retaylor(EKT,TaylorL,Retaylor)       !̩��΢�߶��µ�������ŵ��
!     CALL CAL_Kolmogorov(DISSIP,KolmL)             !����Kolmogorov�߶�
!     CALL CAL_IntegralL(EKT,DISSIP,IntegrL )       !������ֳ߶�
!     CALL CAL_TurnoverTime(EKT,IntegrL ,Turntime)  !��������תʱ��
!����֤���ϳ����������
! ---&---1---------2---------3---------4---------5---------6---------7--
	MODULE FFT
	INTEGER,SAVE :: NX,NY,NZ,NX2,NY2,NZ2,NXP,NYP,NZP,NX2P,NY2P,NZ2P,KMAX,my_id,NPROC,INDEX_NXYZ,INDEX_Serial
	REAL*8,ALLOCATABLE,SAVE :: WAVX(:),WAVY(:),WAVZ(:),trigx(:),trigy(:),trigZ(:),trigx2(:),trigy2(:),trigZ2(:)
	COMPLEX*16,ALLOCATABLE,SAVE :: CWAVX(:),CWAVY(:),CWAVZ(:)
	REAL*8,ALLOCATABLE,SAVE :: EVK2(:,:,:),VECTORK(:,:,:,:)
	END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
	MODULE STATS_SPECTRUM
	REAL*8,SAVE ::  EKT,EKT0,DISSIP,TaylorL,Retaylor,KolmL,IntegrL,Turntime 
	REAL*8,ALLOCATABLE,SAVE :: EKP(:),EK(:),EKX(:),EKXYZ(:,:),EKG(:),KG(:) 
	REAL*8,SAVE :: EKTXYZ(1:3)
	END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
	MODULE STATS_PHYSICAL
	REAL*8,SAVE :: VT2P(1:3),VT2(1:3)  
! 	REAL*8,ALLOCATABLE,SAVE ::                       
	END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
	MODULE VELOCITY
	COMPLEX*16, ALLOCATABLE,SAVE :: V(:,:,:,:),VCW(:,:,:,:),VCWm(:,:,:,:),Wo(:,:,:,:)
	COMPLEX*16,SAVE :: VCWDOTK
	REAL*8, ALLOCATABLE,SAVE :: Vp(:,:,:,:)       !����ռ��ٶ�
	END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
	PROGRAM ISODNS3D

	USE FFT
	USE STATS_SPECTRUM
	USE STATS_PHYSICAL
    USE VELOCITY
	INCLUDE 'com_iso_3d.h'

include 'mpif.h'
integer status(MPI_Status_Size)
INTEGER ierr

	character*2,external :: itochar
    REAL*8 K0,A

	REAL*8  b(0:1)
	INTEGER NT
	REAL*8  TIME,TIME1,TIME2,TIME3
	INTEGER I,J,K,L,it,kx,ky,kz

    REAL*8 FORCE                             !����ϵ��

	REAL*8, ALLOCATABLE :: FRX(:)            !X�����������غ���FX(r)
	REAL*8, ALLOCATABLE :: GRX(:)            !X������������غ���GX(r)
	REAL*8, ALLOCATABLE :: HRX(:)            !X������������غ���HX(r)
	REAL*8, ALLOCATABLE :: GFR(:)            !X������������غ���GF(r)
	REAL*8, ALLOCATABLE :: FRY(:)            !Y������������غ���FY(r)
	REAL*8, ALLOCATABLE :: GRY(:)            !Y������������غ���GY(r)
	REAL*8, ALLOCATABLE :: HRY(:)            !Y������������غ���HY(r)
	REAL*8, ALLOCATABLE :: FRZ(:)            !Z������������غ���FZ(r)
	REAL*8, ALLOCATABLE :: GRZ(:)            !Z������������غ���GZ(r)
	REAL*8, ALLOCATABLE :: HRZ(:)            !Z������������غ���HZ(r)

	REAL*8  SKE(1:4),FLAT(1:4)               !����ռ��ٶȳ���ƫ�Ⱥ�ƽ̹��
    REAL*8  SKEDV(1:4),FLATDV(1:4)           !�ٶ������ʵ�ƫ�Ⱥ�ƽ̹��          
   
! 	COMPLEX*16,ALLOCATABLE :: Div_V(:,:,:)



call MPI_Init(ierr)
call MPI_Comm_Rank(MPI_Comm_World,my_id,ierr)
call MPI_Comm_Size(MPI_Comm_World,NPROC,ierr)
! 	my_id = 0;	NPROC = 1

! ---
IF(my_id.EQ.0) THEN
	OPEN(10,FILE='control_iso_3d.dat',STATUS="UNKNOWN")
	    READ(10,*)
		READ(10,*) NX,NY,NZ
		WRITE(*,*) NX,NY,NZ
		READ(10,*)
		READ(10,*) Re,DT
		WRITE(*,*) Re,DT
		READ(10,*) 
		READ(10,*) K0,A
		WRITE(*,*) K0,A
		READ(10,*)
		READ(10,*) XL,YL,ZL,NSTEPS,IDUMFRQ,IPRNFRQ,IWRITFRQ
		WRITE(*,*) XL,YL,ZL,NSTEPS,IDUMFRQ,IPRNFRQ,IWRITFRQ
		READ(10,*)
		READ(10,*) ICONT,IFORCE,KFMAX,INDEX_NXYZ,INDEX_Serial
		WRITE(*,*) ICONT,IFORCE,KFMAX,INDEX_NXYZ,INDEX_Serial
	CLOSE(10)
ENDIF
                      
CALL MPI_BCAST(NX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(NY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(NZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(Re,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(DT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(K0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(A,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(XL,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(YL,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ZL,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

CALL MPI_BCAST(NSTEPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(IDUMFRQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(IPRNFRQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(IWRITFRQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

CALL MPI_BCAST(ICONT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(IFORCE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(KFMAX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(INDEX_NXYZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(INDEX_Serial,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


	NX2 = 3*NX/2;		NY2 = 3*NY/2;		NZ2 = 3*NZ/2
	NXP = NX/NPROC;		NYP = NY/NPROC;		NZP = NZ/NPROC
	NX2P = NX2/NPROC;	NY2P = NY2/NPROC;	NZ2P = NZ2/NPROC
	
    KMAX=NX/3
	ALLOCATE (Vp(3,0:NX-1,0:NYP-1,0:NZ-1))
	ALLOCATE(FRX(0:NX/2),GRX(0:NX/2),HRX(0:NX/2),GFR(0:NX/2),FRY(0:NY/2),GRY(0:NY/2),HRY(0:NY/2),FRZ(0:NZ/2),GRZ(0:NZ/2),HRZ(0:NZ/2))

    ALLOCATE(EKP(0:KMAX),EK(0:KMAX),EKX(1:KMAX),KG(1:KMAX),EKG(1:KMAX))

! 	ALLOCATE (Div_V(0:NX/2,0:NY-1,0:NZP-1))

	ALLOCATE (V(3,0:NX/2,0:NY-1,0:NZP-1),VCW(3,0:NX/2,0:NY-1,0:NZP-1),VCWm(3,0:NX/2,0:NY-1,0:NZP-1),Wo(1:3,0:NX/2,0:NY-1,0:NZP-1))
	ALLOCATE (WAVX(0:NX/2),WAVY(0:NY-1),WAVZ(0:NZP-1),trigx(2*NX+15),trigy(4*NY+15), &
					trigz(4*NZ+15),trigx2(2*NX2+15),trigy2(4*NY2+15),trigz2(4*NZ2+15) )

	ALLOCATE (CWAVX(0:NX/2),CWAVY(0:NY-1),CWAVZ(0:NZP-1))

	ALLOCATE (EVK2(0:NX/2,0:NY-1,0:NZP-1),VECTORK(3,0:NX/2,0:NY-1,0:NZP-1))



	CALL CPU_TIME (time1)
	CALL SETUP 
	IF(ICONT==0) THEN
		CALL INITIAL (V,A,K0)
		NT = 0;		TIME = 0.0D0
	ELSE
		OPEN(51,FILE='field='//itochar(my_id)//'.dat',STATUS='UNKNOWN')
		READ(51,*) NT,TIME
  		READ(51,*) Vp
		CLOSE(51)
		DO I = 1,3
			CALL xyzfft ('F',Vp(I,0:NX-1,0:NYP-1,0:NZ-1),V(I,0:NX/2,0:NY-1,0:NZP-1))
		ENDDO
	ENDIF
	
    IF(my_id.eq.0) THEN
       CALL OPENFIELD  !CHANGE ON 5.13
	ENDIF



!��ʼʱ���ƽ�����
 DO it = 0,NSTEPS
	IF(it.NE.0) THEN  
		V(1:3,0,0,0)=0.0D0           !��֤����ռ�ƽ���ٶȵ���0

		IF(it<3) THEN
			b(0) = DFLOAT(1+it)/2.0D0;  b(1) = DFLOAT(1-it)/2.0D0
			b = DT*b
		ENDIF
	  
		TIME = TIME + DT

		CALL Curl (Wo,V)
		CALL CROSSTIME2 (VCW,V,Wo)  
		DO kx = 0,NX/2
			DO ky = 0,NY-1
				DO kz = 0,NZP-1
					VCWDOTK = VCW(1,kx,ky,kz)*VECTORK(1,kx,ky,kz) + &
							+ VCW(2,kx,ky,kz)*VECTORK(2,kx,ky,kz) + &
							+ VCW(3,kx,ky,kz)*VECTORK(3,kx,ky,kz)
					VCW(1:3,kx,ky,kz) = VCW(1:3,kx,ky,kz) - VCWDOTK*VECTORK(1:3,kx,ky,kz)
					V(1:3,kx,ky,kz) = ( V(1:3,kx,ky,kz) + ( b(0)*VCW (1:3,kx,ky,kz) + &
															b(1)*VCWm(1:3,kx,ky,kz)/EVK2(kx,ky,kz) ) )/EVK2(kx,ky,kz)
				ENDDO
			ENDDO
		ENDDO
		VCWm = VCW

		IF(IFORCE.EQ.1) THEN
		   CALL CAL_FORCE(V)
		ENDIF

                                                 !ʲôҲ��������������ʼ�ٶȳ��ĸ���ͳ������
	ENDIF

!�����׿ռ�ͳ����
	CALL CAL_SPECTRUM

	DO I=1,3
	   CALL xyzfft ('B',Vp(I,0:NX-1,0:NYP-1,0:NZ-1),V(I,0:NX/2,0:NY-1,0:NZP-1))
	ENDDO
	CALL CALP_ENERGY(VP) 

	    IF(my_id.EQ.0) THEN
	        
			WRITE(100,'(8(F10.5))')TIME,EKT,DISSIP,TaylorL,Retaylor,KolmL,IntegrL,Turntime
        	
			IF(MOD(it,100)==0) THEN
	           WRITE(800,*) 'zone T="time',TIME,'"'
	           DO I=1,KMAX
	              WRITE(800,*) I,EK(I)
	           ENDDO

	           WRITE(900,*) 'zone T="time',TIME,'"'
	           DO I=1,KMAX
	              WRITE(900,*) I,EKX(I)
	           ENDDO

		       WRITE(1000,*)'zone T="time',TIME,'"'
	           DO I=1,KMAX
	              WRITE(1000,*)KG(I),EKG(I)
               ENDDO
            ENDIF 
			WRITE(1100,'(8(F10.5))')TIME,VT2(1),VT2(2),VT2(3),0.5d0*(VT2(1)+VT2(2)+VT2(3))
		    
			WRITE(*,*) 'it = ', it,'    TIME = ', TIME
		    WRITE(*,*)"EKT,DISSIP,TaylorL,Kolml,IntegrL"
        	WRITE(*,*)EKT,DISSIP,TaylorL,Kolml,IntegrL
	        WRITE(*,*)"Retaylor,Turntime"
	        WRITE(*,*)Retaylor,Turntime
	        write(*,*)'V(1,0,0,0)=',V(1,0,0,0)
	        write(*,*)'V(2,0,0,0)=',V(2,0,0,0)
        	WRITE(*,*)'EK(1)=',EK(1)
	        WRITE(*,*)'EK(2)=',EK(2)
        	WRITE(*,*)'EK(KMAX-1)=',EK(KMAX-1)
	        WRITE(*,*)'EK(KMAX)=',EK(KMAX-1)
            WRITE(*,*)"-----------------------"
	    ENDIF  
	    
		IF(MOD(it+1,200)==0) THEN
			DO I=1,3
				CALL xyzfft ('B',Vp(I,0:NX-1,0:NYP-1,0:NZ-1),V(I,0:NX/2,0:NY-1,0:NZP-1))
			ENDDO
            
			OPEN(61,FILE='field='//itochar(my_id)//'.dat',STATUS='UNKNOWN')       !'//itochar(my_id)//'
 			WRITE(61,*) NT+it,TIME
  			WRITE(61,*)Vp(1:3,0:NX-1,0:NYP-1,0:NZ-1)
			CLOSE(61)
		ENDIF






ENDDO   !������ѭ��
 	CALL CPU_TIME (time2)

	CALL MPI_FINALIZE(ierr)
	WRITE(*,*) 'Time = ', time2-time1     
	END



SUBROUTINE OPENFIELD
!����ռ�ͳ������غ���F(r),G(r)
!���д�ļ����кܴ�����
!-------------------------------------------------------
! 	        OPEN (20,FILE='Fr.plt',status='unknown')  
! 		    WRITE(20,*)'VARIABLES="r","frx","gry","hrz"'
! 		    OPEN (30,FILE='Gr.plt',status='unknown')
! 		    WRITE(30,*)'VARIABLES="r","grx","gfr","fry","hry","frz","grz","hrx"'
! 			OPEN (40,FILE='Skewness.plt',status='unknown')
! 			WRITE(40,*)'VARIABLES="time","skewnessU","skewnessV","skewnessW","skewnessUmean"'
! 			OPEN (50,FILE='flatness.plt',status='unknown')
! 			WRITE(50,*)'VARIABLES="time","flatnessU","flatnessV","flatnessW","flatnessUmean"'	
! 			OPEN (60,FILE='SkewDU.plt',status='unknown')
! 			WRITE(60,*)'VARIABLES="time","skewnessDU","skewnessDV","skewnessDW","skewnessDUmean"'
! 			OPEN (70,FILE='flatDU.plt',status='unknown')
! 			WRITE(70,*)'VARIABLES="time","flatnessDU","flatnessDV","flatnessDW","flatnessDUmean"'
!-------------------------------------------------------
!������ 5.17 18:00


!�׿ռ�ͳ����
			OPEN (100,FILE='total energy.plt',status='unknown')
			WRITE(100,*)'VARIABLES="time","EKT","Dissip","TaylorL","Retaylor","KolmL","L","T"'

	
			OPEN (800,FILE='3Denergy spectrum.plt',status='unknown')
			WRITE(800,*)'VARIABLES="wavenumber","energy spectrum"'
			OPEN (900,FILE='1D X energy spectrum.plt',status='unknown')
			WRITE(900,*)'VARIABLES="wavenumber","energy spectrum"'
			OPEN (1000,FILE='GuiYi energy spcetrum.plt',status='unknown')
			WRITE(1000,*)'VARIABLES="guiyi K","guiyi EK"'

			OPEN (1100,FILE='U2V2W2.plt',status='unknown')
			WRITE(1100,*)'VARIABLES="time","U2","V2","W2","0.5UVW2"'
! 
! 			OPEN (1200,FILE='CHECK X energy spcetrum.plt',status='unknown')
! 			WRITE(1200,*)'VARIABLES="TIME","EKTX"'
! 			OPEN (1300,FILE='CHECK Y energy spcetrum.plt',status='unknown')
! 			WRITE(1300,*)'VARIABLES="TIME","EKTY"'
! 			OPEN (1400,FILE='CHECK Z energy spcetrum.plt',status='unknown')
! 			WRITE(1400,*)'VARIABLES="TIME","EKTZ"'
RETURN
END