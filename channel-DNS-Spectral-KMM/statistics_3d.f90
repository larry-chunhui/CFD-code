! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE statics_av (savq,RrEx,RrEz,Vp,PREp)

	  USE FFT
	  INCLUDE 'dim.h'

	  REAL*8, DIMENSION(0:NX-1,0:NY,0:NZP-1) :: up,vp,wp,prep,Wxp,Wyp,Wzp
	  REAL*8 savq(49,0:NY),RrEx(0:NX/2,35,0:NY),RrEz(0:NZ/2,35,0:NY)
	  INTEGER idx,ix2,idz,iz2

	  REAL*8, ALLOCATABLE :: Exx(:),Ezz(:)

      ALLOCATE(Exx(0:NX/2),Ezz(0:NZ/2))


!	  REAL*8, DIMENSION(49,0:NYP) :: savq

!savq includ u,v,w,p,Wx,Wy,Wz;
!            u^2,v^2,w^2,p^2,Wx^2,Wy^2,Wz^2;
!            u^3,v^3,w^3,p^3,Wx^3,Wy^3,Wz^3;
!            u^4,v^4,w^4,p^4,Wx^4,Wy^4,Wz^4;
!            uv,vw,wu,WxWy,WyWz,WzWx,pu,pv,pw;
!            uWx,uWy,uWz,vWx,vWy,vWz,wWx,wWy,wWz,pWx,pWy,pWz

!	  REAL*8, DIMENSION(0:NX/2,35,0:NYP) :: RrEx
!	  REAL*8, DIMENSION(0:NZ/2,35,0:NYP) :: RrEz

!RrEx includ uux,vvx,wwx,ppx,WxWxx,WyWyx,WzWzx;
!            uvx,vwx,wux,WxWyx,WyWzx,WzWxx,pux,pvx,pwx;
!            uWxx,uWyx,uWzx,vWxx,vWyx,vWzx,wWxx,wWyx,wWzx,pWxx,pWyx,pWzx
!            Euux,Evvx,Ewwx,Eppx,EWxWxx,EWyWyx,EWzWzx;
!RrEz includ uuz,vvz,wwz,ppz,WxWxz,WyWyz,WzWzz;
!            uvz,vwz,wuz,WxWyz,WyWzz,WzWxz,puz,pvz,pwz;
!            uWxz,uWyz,uWzz,vWxz,vWyz,vWzz,wWxz,wWyz,wWzz,pWxz,pWyz,pWzz
!            Euuz,Evvz,Ewwz,Eppz,EWxWxz,EWyWyz,EWzWzz

	  savq = 0.0D0
	  RrEx = 0.0D0
	  RrEz = 0.0D0

      DO iy=0,NY

         DO iz=0,NZ-1
            DO ix=0,NX-1

		       savq( 1,iy) = savq( 1,iy) + up(ix,iy,iz)
		       savq( 2,iy) = savq( 2,iy) + vp(ix,iy,iz)
		       savq( 3,iy) = savq( 3,iy) + wp(ix,iy,iz)
		       savq( 4,iy) = savq( 4,iy) + prep(ix,iy,iz)
		       savq( 5,iy) = savq( 5,iy) + Wxp(ix,iy,iz)
		       savq( 6,iy) = savq( 6,iy) + Wyp(ix,iy,iz)
		       savq( 7,iy) = savq( 7,iy) + Wzp(ix,iy,iz)

		       savq( 8,iy) = savq( 8,iy) + up(ix,iy,iz)**2
		       savq( 9,iy) = savq( 9,iy) + vp(ix,iy,iz)**2
		       savq(10,iy) = savq(10,iy) + wp(ix,iy,iz)**2
		       savq(11,iy) = savq(11,iy) + prep(ix,iy,iz)**2
		       savq(12,iy) = savq(12,iy) + Wxp(ix,iy,iz)**2
		       savq(13,iy) = savq(13,iy) + Wyp(ix,iy,iz)**2
		       savq(14,iy) = savq(14,iy) + Wzp(ix,iy,iz)**2

		       savq(15,iy) = savq(15,iy) + up(ix,iy,iz)**3
		       savq(16,iy) = savq(16,iy) + vp(ix,iy,iz)**3
		       savq(17,iy) = savq(17,iy) + wp(ix,iy,iz)**3
		       savq(18,iy) = savq(18,iy) + prep(ix,iy,iz)**3
		       savq(19,iy) = savq(19,iy) + Wxp(ix,iy,iz)**3
		       savq(20,iy) = savq(20,iy) + Wyp(ix,iy,iz)**3
		       savq(21,iy) = savq(21,iy) + Wzp(ix,iy,iz)**3

		       savq(22,iy) = savq(22,iy) + up(ix,iy,iz)**4
		       savq(23,iy) = savq(23,iy) + vp(ix,iy,iz)**4
		       savq(24,iy) = savq(24,iy) + wp(ix,iy,iz)**4
		       savq(25,iy) = savq(25,iy) + prep(ix,iy,iz)**4
		       savq(26,iy) = savq(26,iy) + Wxp(ix,iy,iz)**4
		       savq(27,iy) = savq(27,iy) + Wyp(ix,iy,iz)**4
		       savq(28,iy) = savq(28,iy) + Wzp(ix,iy,iz)**4

		       savq(29,iy) = savq(29,iy) + up(ix,iy,iz)*vp(ix,iy,iz)
		       savq(30,iy) = savq(30,iy) + vp(ix,iy,iz)*wp(ix,iy,iz)
		       savq(31,iy) = savq(31,iy) + wp(ix,iy,iz)*up(ix,iy,iz)
		       savq(32,iy) = savq(32,iy) + Wxp(ix,iy,iz)*Wyp(ix,iy,iz)
		       savq(33,iy) = savq(33,iy) + Wyp(ix,iy,iz)*Wzp(ix,iy,iz)
		       savq(34,iy) = savq(34,iy) + Wzp(ix,iy,iz)*Wxp(ix,iy,iz)
		       savq(35,iy) = savq(35,iy) + prep(ix,iy,iz)*up(ix,iy,iz)
		       savq(36,iy) = savq(36,iy) + prep(ix,iy,iz)*vp(ix,iy,iz)
		       savq(37,iy) = savq(37,iy) + prep(ix,iy,iz)*wp(ix,iy,iz)

		       savq(38,iy) = savq(38,iy) + up(ix,iy,iz)*Wxp(ix,iy,iz)
		       savq(39,iy) = savq(39,iy) + up(ix,iy,iz)*Wyp(ix,iy,iz)
		       savq(40,iy) = savq(40,iy) + up(ix,iy,iz)*Wzp(ix,iy,iz)
		       savq(41,iy) = savq(41,iy) + vp(ix,iy,iz)*Wxp(ix,iy,iz)
		       savq(42,iy) = savq(42,iy) + vp(ix,iy,iz)*Wyp(ix,iy,iz)
		       savq(43,iy) = savq(43,iy) + vp(ix,iy,iz)*Wzp(ix,iy,iz)
		       savq(44,iy) = savq(44,iy) + wp(ix,iy,iz)*Wxp(ix,iy,iz)
		       savq(45,iy) = savq(45,iy) + wp(ix,iy,iz)*Wyp(ix,iy,iz)
		       savq(46,iy) = savq(46,iy) + wp(ix,iy,iz)*Wzp(ix,iy,iz)
		       savq(47,iy) = savq(47,iy) + prep(ix,iy,iz)*Wxp(ix,iy,iz)
		       savq(48,iy) = savq(48,iy) + prep(ix,iy,iz)*Wyp(ix,iy,iz)
		       savq(49,iy) = savq(49,iy) + prep(ix,iy,iz)*Wzp(ix,iy,iz)
	     
               DO idx=0,NX/2
		          ix2 = mod(ix+idx,NX)

                  RrEx(idx, 1,iy) = RrEx(idx, 1,iy) + up(ix,iy,iz)*up(ix2,iy,iz)
                  RrEx(idx, 2,iy) = RrEx(idx, 2,iy) + vp(ix,iy,iz)*vp(ix2,iy,iz)
                  RrEx(idx, 3,iy) = RrEx(idx, 3,iy) + wp(ix,iy,iz)*wp(ix2,iy,iz)
                  RrEx(idx, 4,iy) = RrEx(idx, 4,iy) + prep(ix,iy,iz)*prep(ix2,iy,iz)
                  RrEx(idx, 5,iy) = RrEx(idx, 5,iy) + Wxp(ix,iy,iz)*Wxp(ix2,iy,iz)
                  RrEx(idx, 6,iy) = RrEx(idx, 6,iy) + Wyp(ix,iy,iz)*Wyp(ix2,iy,iz)
                  RrEx(idx, 7,iy) = RrEx(idx, 7,iy) + Wzp(ix,iy,iz)*Wzp(ix2,iy,iz)

                  RrEx(idx, 8,iy) = RrEx(idx, 8,iy) + up(ix,iy,iz)*vp(ix2,iy,iz)
                  RrEx(idx, 9,iy) = RrEx(idx, 9,iy) + vp(ix,iy,iz)*wp(ix2,iy,iz)
                  RrEx(idx,10,iy) = RrEx(idx,10,iy) + wp(ix,iy,iz)*up(ix2,iy,iz)
                  RrEx(idx,11,iy) = RrEx(idx,11,iy) + Wxp(ix,iy,iz)*Wyp(ix2,iy,iz)
                  RrEx(idx,12,iy) = RrEx(idx,12,iy) + Wyp(ix,iy,iz)*Wzp(ix2,iy,iz)
                  RrEx(idx,13,iy) = RrEx(idx,13,iy) + Wzp(ix,iy,iz)*Wxp(ix2,iy,iz)
                  RrEx(idx,14,iy) = RrEx(idx,14,iy) + prep(ix,iy,iz)*up(ix2,iy,iz)
                  RrEx(idx,15,iy) = RrEx(idx,15,iy) + prep(ix,iy,iz)*vp(ix2,iy,iz)
                  RrEx(idx,16,iy) = RrEx(idx,16,iy) + prep(ix,iy,iz)*wp(ix2,iy,iz)

                  RrEx(idx,17,iy) = RrEx(idx,17,iy) + up(ix,iy,iz)*Wxp(ix2,iy,iz)
                  RrEx(idx,18,iy) = RrEx(idx,18,iy) + up(ix,iy,iz)*Wyp(ix2,iy,iz)
                  RrEx(idx,19,iy) = RrEx(idx,19,iy) + up(ix,iy,iz)*Wzp(ix2,iy,iz)
                  RrEx(idx,20,iy) = RrEx(idx,20,iy) + vp(ix,iy,iz)*Wxp(ix2,iy,iz)
                  RrEx(idx,21,iy) = RrEx(idx,21,iy) + vp(ix,iy,iz)*Wyp(ix2,iy,iz)
                  RrEx(idx,22,iy) = RrEx(idx,22,iy) + vp(ix,iy,iz)*Wzp(ix2,iy,iz)
                  RrEx(idx,23,iy) = RrEx(idx,23,iy) + wp(ix,iy,iz)*Wxp(ix2,iy,iz)
                  RrEx(idx,24,iy) = RrEx(idx,24,iy) + wp(ix,iy,iz)*Wyp(ix2,iy,iz)
                  RrEx(idx,25,iy) = RrEx(idx,25,iy) + wp(ix,iy,iz)*Wzp(ix2,iy,iz)
                  RrEx(idx,26,iy) = RrEx(idx,26,iy) + prep(ix,iy,iz)*Wxp(ix2,iy,iz)
                  RrEx(idx,27,iy) = RrEx(idx,27,iy) + prep(ix,iy,iz)*Wyp(ix2,iy,iz)
                  RrEx(idx,28,iy) = RrEx(idx,28,iy) + prep(ix,iy,iz)*Wzp(ix2,iy,iz)

               ENDDO

               DO idz=0,NZ/2
		          iz2 = mod(iz+idz,NZ)

                  RrEz(idz, 1,iy) = RrEz(idz, 1,iy) + up(ix,iy,iz)*up(ix,iy,iz2)
                  RrEz(idz, 2,iy) = RrEz(idz, 2,iy) + vp(ix,iy,iz)*vp(ix,iy,iz2)
                  RrEz(idz, 3,iy) = RrEz(idz, 3,iy) + wp(ix,iy,iz)*wp(ix,iy,iz2)
                  RrEz(idz, 4,iy) = RrEz(idz, 4,iy) + prep(ix,iy,iz)*prep(ix,iy,iz2)
                  RrEz(idz, 5,iy) = RrEz(idz, 5,iy) + Wxp(ix,iy,iz)*Wxp(ix,iy,iz2)
                  RrEz(idz, 6,iy) = RrEz(idz, 6,iy) + Wyp(ix,iy,iz)*Wyp(ix,iy,iz2)
                  RrEz(idz, 7,iy) = RrEz(idz, 7,iy) + Wzp(ix,iy,iz)*Wzp(ix,iy,iz2)

                  RrEz(idz, 8,iy) = RrEz(idz, 8,iy) + up(ix,iy,iz)*vp(ix,iy,iz2)
                  RrEz(idz, 9,iy) = RrEz(idz, 9,iy) + vp(ix,iy,iz)*wp(ix,iy,iz2)
                  RrEz(idz,10,iy) = RrEz(idz,10,iy) + wp(ix,iy,iz)*up(ix,iy,iz2)
                  RrEz(idz,11,iy) = RrEz(idz,11,iy) + Wxp(ix,iy,iz)*Wyp(ix,iy,iz2)
                  RrEz(idz,12,iy) = RrEz(idz,12,iy) + Wyp(ix,iy,iz)*Wzp(ix,iy,iz2)
                  RrEz(idz,13,iy) = RrEz(idz,13,iy) + Wzp(ix,iy,iz)*Wxp(ix,iy,iz2)
                  RrEz(idz,14,iy) = RrEz(idz,14,iy) + prep(ix,iy,iz)*up(ix,iy,iz2)
                  RrEz(idz,15,iy) = RrEz(idz,15,iy) + prep(ix,iy,iz)*vp(ix,iy,iz2)
                  RrEz(idz,16,iy) = RrEz(idz,16,iy) + prep(ix,iy,iz)*wp(ix,iy,iz2)

                  RrEz(idz,17,iy) = RrEz(idz,17,iy) + up(ix,iy,iz)*Wxp(ix,iy,iz2)
                  RrEz(idz,18,iy) = RrEz(idz,18,iy) + up(ix,iy,iz)*Wyp(ix,iy,iz2)
                  RrEz(idz,19,iy) = RrEz(idz,19,iy) + up(ix,iy,iz)*Wzp(ix,iy,iz2)
                  RrEz(idz,20,iy) = RrEz(idz,20,iy) + vp(ix,iy,iz)*Wxp(ix,iy,iz2)
                  RrEz(idz,21,iy) = RrEz(idz,21,iy) + vp(ix,iy,iz)*Wyp(ix,iy,iz2)
                  RrEz(idz,22,iy) = RrEz(idz,22,iy) + vp(ix,iy,iz)*Wzp(ix,iy,iz2)
                  RrEz(idz,23,iy) = RrEz(idz,23,iy) + wp(ix,iy,iz)*Wxp(ix,iy,iz2)
                  RrEz(idz,24,iy) = RrEz(idz,24,iy) + wp(ix,iy,iz)*Wyp(ix,iy,iz2)
                  RrEz(idz,25,iy) = RrEz(idz,25,iy) + wp(ix,iy,iz)*Wzp(ix,iy,iz2)
                  RrEz(idz,26,iy) = RrEz(idz,26,iy) + prep(ix,iy,iz)*Wxp(ix,iy,iz2)
                  RrEz(idz,27,iy) = RrEz(idz,27,iy) + prep(ix,iy,iz)*Wyp(ix,iy,iz2)
                  RrEz(idz,28,iy) = RrEz(idz,28,iy) + prep(ix,iy,iz)*Wzp(ix,iy,iz2)

              ENDDO

           ENDDO
		 ENDDO

		 call Exz (NX,NZ,Exx,Ezz,up(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,29,iy) = RrEx(0:NX/2,29,iy) + Exx
         RrEz(0:NZ/2,29,iy) = RrEz(0:NZ/2,29,iy) + Ezz
		 call Exz (NX,NZ,Exx,Ezz,vp(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,30,iy) = RrEx(0:NX/2,30,iy) + Exx
         RrEz(0:NZ/2,30,iy) = RrEz(0:NZ/2,30,iy) + Ezz
		 call Exz (NX,NZ,Exx,Ezz,wp(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,31,iy) = RrEx(0:NX/2,31,iy) + Exx
         RrEz(0:NZ/2,31,iy) = RrEz(0:NZ/2,31,iy) + Ezz
		 call Exz (NX,NZ,Exx,Ezz,prep(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,32,iy) = RrEx(0:NX/2,32,iy) + Exx
         RrEz(0:NZ/2,32,iy) = RrEz(0:NZ/2,32,iy) + Ezz
		 call Exz (NX,NZ,Exx,Ezz,Wxp(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,33,iy) = RrEx(0:NX/2,33,iy) + Exx
         RrEz(0:NZ/2,33,iy) = RrEz(0:NZ/2,33,iy) + Ezz
		 call Exz (NX,NZ,Exx,Ezz,Wyp(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,34,iy) = RrEx(0:NX/2,34,iy) + Exx
         RrEz(0:NZ/2,34,iy) = RrEz(0:NZ/2,34,iy) + Ezz
		 call Exz (NX,NZ,Exx,Ezz,Wzp(0:NX-1,iy,0:NZ-1),trigx,trigzr)
         RrEx(0:NX/2,35,iy) = RrEx(0:NX/2,35,iy) + Exx
         RrEz(0:NZ/2,35,iy) = RrEz(0:NZ/2,35,iy) + Ezz

	  ENDDO
	  savq = savq/DFLOAT(NZ*NX)
	  RrEx = RrEx/DFLOAT(NZ*NX)
	  RrEz = RrEz/DFLOAT(NZ*NX)

	  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE Exz (NX,NZ,Exx,Ezz,UM,trigx,trigzr)

	  IMPLICIT REAL*8 (A-H, O-Z)

	  REAL*8 UM(0:NX-1,0:NZ-1),Exx(0:NX/2),Ezz(0:NZ/2)

	  REAL*8, ALLOCATABLE :: UMX(:),UMZ(:)
	  COMPLEX*16, ALLOCATABLE :: CEx(:),CEz(:)

      ALLOCATE(UMX(0:NX+1),UMZ(0:NZ+1),CEx(0:NX/2),CEz(0:NZ/2))

	  Exx = 0.0D0
	  DO iz = 0,NZ-1
	     UMX(0:NX-1) = UM(0:NX-1,iz)
		 CALL rfftf(NX,UMX,trigx)
	     CEx(0) = DCMPLX(UMX(0),0.0D0)
		 DO i=1,NX/2-1
	        CEx(i) = DCMPLX(UMX(2*i-1),UMX(2*i))
		 ENDDO
	     CEx(NX/2) = DCMPLX(UMX(NX-1),0.0D0)
	     Exx = Exx + DREAL(CEx*DCONJG(CEx))/DFLOAT(NX)
	  ENDDO

	  Ezz = 0.0D0
	  DO ix = 0,NX-1
	     UMZ(0:NZ-1) = UM(ix,0:NZ-1)
		 CALL rfftf(NZ,UMZ,trigz0)
	     CEz(0) = DCMPLX(UMZ(0),0.0D0)
		 DO i=1,NZ/2-1
	        CEz(i) = DCMPLX(UMZ(2*i-1),UMZ(2*i))
		 ENDDO
	     CEz(NZ/2) = DCMPLX(UMZ(NZ-1),0.0D0)
	     Ezz = Ezz + DREAL(CEz*DCONJG(CEz))/DFLOAT(NZ)
	  ENDDO

	  END