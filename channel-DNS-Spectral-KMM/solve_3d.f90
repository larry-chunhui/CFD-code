! ---&---1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE SOLVE (U,F,G,a,b,BC,WN,AD)

	  USE FFT

	  IMPLICIT NONE

      COMPLEX*16, DIMENSION(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1) :: U,F
      COMPLEX*16, DIMENSION(0:NXP/2-1,0:1,-NZ/2+1:NZ/2-1) :: BC

      REAL*8  WN(0:NXP/2-1,-NZ/2+1:NZ/2-1),AD(0:NXP/2-1,0:NY,-NZ/2+1:NZ/2-1)

      COMPLEX*16 COR

	  REAL*8 G,a,b,FUNC

	  INTEGER iy,index,IBF

	  IF(G==0.0D0 .AND. a==0.0D0 .AND. my_id==0) THEN
	     IF(b==0.0D0) THEN
			COR=0.0D0
		 ELSE
			COR = (BC(0,0,0)+BC(0,1,0))/2.0D0/b
		 ENDIF
	     DO iy=2,NY,2
	        COR = COR + F(0,iy,0)/DFLOAT(iy**2-1)
	     ENDDO

! 	     F(0,0,0) = COR

		 COR = COR - F(0,0,0)
	     BC(0,0:1,0) = BC(0,0:1,0) - b*COR
	  ENDIF

	  
	  DO iy = 2,NY-2    !æÿ’Û”“∂ÀœÓ
         U(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = WL(iy)*F(0:NXP/2-1,iy-2,-NZ/2+1:NZ/2-1) &
										+ WD(iy)*F(0:NXP/2-1, iy ,-NZ/2+1:NZ/2-1) &
										+ WR(iy)*F(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1)
      ENDDO
      U(0:NXP/2-1,ny-1,-NZ/2+1:NZ/2-1) = WL(ny-1)*F(0:NXP/2-1,ny-3,-NZ/2+1:NZ/2-1) &
										+ WD(ny-1)*F(0:NXP/2-1,ny-1,-NZ/2+1:NZ/2-1)
      U(0:NXP/2-1, ny ,-NZ/2+1:NZ/2-1) = WL( ny )*F(0:NXP/2-1,ny-2,-NZ/2+1:NZ/2-1) &
										+ WD( ny )*F(0:NXP/2-1, ny ,-NZ/2+1:NZ/2-1)

      U(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)=(BC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)+BC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1) )/2.0D0
      U(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)=(BC(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)-BC(0:NXP/2-1,1,-NZ/2+1:NZ/2-1) )/2.0D0

      DO iy = NY-2,2,-1
         U(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = U(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) - WR(iy)*WN &
									*AD(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1)*U(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1)
      ENDDO

       DO iy = NY,NY-1,-1
         index = MOD(iy,2)
		 IF(a==0.0D0 .AND. b==0.0D0) THEN
			IBF=1-INDEX
			FUNC = DFLOAT(IBF*iy**2+1-IBF)
		 ELSE
			FUNC = a+b*DFLOAT(iy**2)
		 ENDIF
         AD(0:NXP/2-1,index,-NZ/2+1:NZ/2-1) = FUNC
      ENDDO

      DO iy = NY-2,0,-1
         index = MOD(iy,2)
		 IF(a==0.0D0 .AND. b==0.0D0) THEN
			IBF=1-INDEX
			FUNC = DFLOAT(IBF*iy**2+1-IBF)
		 ELSE
			FUNC = a+b*DFLOAT(iy**2)
		 ENDIF
         U(0:NXP/2-1,index,-NZ/2+1:NZ/2-1)=U(0:NXP/2-1,index,-NZ/2+1:NZ/2-1)-AD(0:NXP/2-1,index,-NZ/2+1:NZ/2-1) &
											*AD(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1)*U(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1)
         AD(0:NXP/2-1,index,-NZ/2+1:NZ/2-1) = FUNC - WL(iy+2)*WN &
										*AD(0:NXP/2-1,index,-NZ/2+1:NZ/2-1)*AD(0:NXP/2-1,iy+2,-NZ/2+1:NZ/2-1)
      ENDDO
	   
	  IF(G==0.0D0 .AND. a==0.0D0 .AND. my_id==0) THEN
	     U(0,0,0) = DCMPLX(0.0D0,0.0D0)
		 AD(0,0,0) = 1.0D0
	  ENDIF

      U(0:NXP/2-1,0,-NZ/2+1:NZ/2-1) = U(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)/AD(0:NXP/2-1,0,-NZ/2+1:NZ/2-1)
      U(0:NXP/2-1,1,-NZ/2+1:NZ/2-1) = U(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)/AD(0:NXP/2-1,1,-NZ/2+1:NZ/2-1)

      DO iy = 2,NY
         U(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) = (U(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1) - WL(iy)*WN &
								*U(0:NXP/2-1,iy-2,-NZ/2+1:NZ/2-1))*AD(0:NXP/2-1,iy,-NZ/2+1:NZ/2-1)
	  ENDDO

      RETURN
      END