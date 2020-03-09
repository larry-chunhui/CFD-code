! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE DIV_c (uvwpc,DVC,G,WNG10,ADG10)

	  USE FFT
	  INCLUDE 'dim.h'

      REAL*8 uvwpc(4,2,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),DVC(2,0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1)
      REAL*8 WNG10(0:NXP/2-1,-NZ/2+1:NZ/2-1),ADG10(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)

	  COMPLEX*16, ALLOCATABLE :: DPREc(:,:,:,:),Vc(:,:,:,:),BC(:,:,:)
	  REAL*8, ALLOCATABLE :: TEMP(:,:,:)

	  ALLOCATE(DPREc(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),Vc(3,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
						BC(0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1),TEMP(2,0:NXP/2-1,-NZ/2+1:NZ/2-1))

	  BC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1) = DCMPLX(1.0D0, 0.0D0)
	  BC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1) = DCMPLX(0.0D0, 1.0D0)
	  DPREc(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) = DCMPLX(0.0D0, 0.0D0)
	  CALL SOLVE_i (WNG10,ADG10,0.0D0,1.0D0,0.0D0)
	  CALL SOLVE (Vc(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),DPREc(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),0.0D0,1.0D0,0.0D0,BC,WNG10,ADG10)

	  DO ix = 0,NXP/2-1
		 DO iy = 0,NY
			DO iz = -NZ/2+1,NZ/2-1
			   uvwpc(4,1,ix,iy,iz) =  DREAL(Vc(1,ix,iy,iz))
			   uvwpc(4,2,ix,iy,iz) =  DIMAG(Vc(1,ix,iy,iz))
			ENDDO
		 ENDDO
	  ENDDO

	  DO I=1,3
		 CALL DERIV (I,DPREc(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),Vc(1,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1))
	  ENDDO

	  DPREc = Re*DPREc

	  BC=DCMPLX(0.0D0,0.0D0)
	  CALL SOLVE_i (WNG10,ADG10,G,1.0D0,0.0D0)
	  DO I=1,3
		 CALL SOLVE (Vc(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1),DPREc(I,0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1), &
											G,1.0D0,0.0D0,BC,WNG10,ADG10)
	  ENDDO

	  DO ix = 0,NXP/2-1
		 DO iy = 0,NY
			DO iz = -NZ/2+1,NZ/2-1
			   uvwpc(1,1,ix,iy,iz) =  DIMAG(Vc(1,ix,iy,iz))
			   uvwpc(1,2,ix,iy,iz) = -DREAL(Vc(1,ix,iy,iz))
			   uvwpc(2,1,ix,iy,iz) =  DREAL(Vc(2,ix,iy,iz))
			   uvwpc(2,2,ix,iy,iz) =  DIMAG(Vc(2,ix,iy,iz))
			   uvwpc(3,1,ix,iy,iz) =  DIMAG(Vc(3,ix,iy,iz))
			   uvwpc(3,2,ix,iy,iz) = -DREAL(Vc(3,ix,iy,iz))
			ENDDO
		 ENDDO
	  ENDDO

	  DVC = 0.0D0
	  DO iy=0,NY
		 TEMP=DFLOAT(iy**2)*uvwpc(2,1:2,0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
		 DVC(1:2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1) = DVC(1:2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1) + TEMP
		 IF(MOD(iy,2)==0) THEN
			DVC(1:2,0:NXP/2-1,1,-NZ/2+1:NZ/2-1) = DVC(1:2,0:NXP/2-1,1,-NZ/2+1:NZ/2-1) - TEMP
		 ELSE
			DVC(1:2,0:NXP/2-1,1,-NZ/2+1:NZ/2-1) = DVC(1:2,0:NXP/2-1,1,-NZ/2+1:NZ/2-1) + TEMP
		 ENDIF
	  ENDDO

	  TEMP(1,0:NXP/2-1,-NZ/2+1:NZ/2-1)=(DVC(1,0:NXP/2-1,0,-NZ/2+1:NZ/2-1)*DVC(2,0:NXP/2-1,1,-NZ/2+1:NZ/2-1) - &
										DVC(1,0:NXP/2-1,1,-NZ/2+1:NZ/2-1)*DVC(2,0:NXP/2-1,0,-NZ/2+1:NZ/2-1))
	  IF(my_id==0) TEMP(1,0,0) = 0.25D0/DVC(1,0,0,0)
  
      DO J=0,1
	     DO I=1,2
		    DVC(I,0:NXP/2-1,J,-NZ/2+1:NZ/2-1)=DVC(I,0:NXP/2-1,J,-NZ/2+1:NZ/2-1)/TEMP(1,0:NXP/2-1,-NZ/2+1:NZ/2-1)
			IF(my_id==0) DVC(I,0,J,0)=TEMP(1,0,0)
		 ENDDO
	  ENDDO

	  DEALLOCATE(DPREc,Vc,BC,TEMP)

	  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE Q_c (UG1,QG1,G)

	  USE FFT
	  INCLUDE 'dim.h'

      REAL*8 UG1(0:NY)

	  REAL*8 WN,FC
	  INTEGER index

	  REAL*8, ALLOCATABLE :: F(:),AD(:)

      ALLOCATE(F(0:NY),AD(0:NY))

	  WN = -G

      AD(ny) = 1.0D0/(1.0D0+ WD(ny)*WN)
      AD(ny-1) = 1.0D0/(1.0D0+ WD(ny-1)*WN)

      DO iy = NY-2,2,-1
         FC = WR(iy)*WL(iy+2)
         AD(iy) = 1.0D0/(1.0D0+(WD(iy)-FC*WN*AD(iy+2))*WN)
      ENDDO

	  F = 0.0D0
	  F(0) = Re
	  
      UG1(0) = 0.0D0
      UG1(1) = 0.0D0
	  DO iy = 2,NY-2
         UG1(iy) = WL(iy)*F(iy-2) &
				 + WD(iy)*F(iy  ) &
				 + WR(iy)*F(iy+2)
      ENDDO
      UG1(ny-1) = WL(ny-1)*F(ny-3) &
				+ WD(ny-1)*F(ny-1)
      UG1( ny ) = WL( ny )*F(ny-2) &
				+ WD( ny )*F( ny )


      DO iy = NY-2,2,-1
         UG1(iy) = UG1(iy) - WR(iy)*WN*AD(iy+2)*UG1(iy+2)
      ENDDO

      DO iy = NY,NY-1,-1
         index = MOD(iy,2)
		 AD(index) = 1.0D0
      ENDDO

      DO iy = NY-2,0,-1
         index = MOD(iy,2)
         UG1(index) = UG1(index) - AD(index)*AD(iy+2)*UG1(iy+2)
         AD(index) = 1.0D0 - WL(iy+2)*WN*AD(index)*AD(iy+2)
      ENDDO	   

      UG1(0) = UG1(0)/AD(0)
      UG1(1) = UG1(1)/AD(1)

      DO iy = 2,NY
         UG1(iy) = (UG1(iy) - WL(iy)*WN*UG1(iy-2))*AD(iy)
	  ENDDO

	  QG1=0.0D0
	  DO iy=0,NY,2
	     QG1 = QG1 - UG1(iy)/DFLOAT(iy**2-1)
	  ENDDO

	  DEALLOCATE(F,AD)

	  END
