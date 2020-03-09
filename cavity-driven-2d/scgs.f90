!c Vanka (1986), JCP, 1986, Vol. 65, p. 138
!c URFU,URFV: ! under-relaxation factor for Vanka matrix solve
program main
USE DFPORT
include 'incl.h'

integer i,j
real*8 time1
! define # of grid points
! define convergence criterion
      MAXIT=100000
      URFV=0.5
      URFU=0.5
      SORMAX=0.0001
! input total number of iterations and Reynolds number
      write(*,*)'MAXIT,RE=?'
      write(*,*)'URFU?'
time1=timef()     
! mesh generation for a 1x1 square cavity
      CALL GRID
! specify initial condition
      CALL INIT
! Vanka's SCGS algorithm
!here is the problem
      CALL SCGS
! post-processing
! calculate U,V at the centroid of P-CV for plotting purpose
	  call output

		time1=timef()
		open(1,file='CPU_time.dat')
		write(1,*) 'CPU_time=',time1
		close(1)
      STOP
      END

!c-------------------------------------------------------
!c-------------------------------------------------------
      SUBROUTINE GRID
      include 'incl.h'
	  integer I,J
      DX=1./FLOAT(NI)  !uniform grid in x and y
      DY=1./FLOAT(NJ)
!c including halo data region
!c
! specify number of "halo-data" layers

      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      X(I,J)=DX*I
      Y(I,J)=DY*J
      END DO
      END DO
      RETURN
      END

!-------------------------------------------------------
! zero velocity && pressure field
!-------------------------------------------------------
      SUBROUTINE INIT
      include 'incl.h'
	  integer I,J
      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      U(I,J)=0.0
      V(I,J)=0.0
      P(I,J)=0.0
      END DO
      END DO
      RETURN
      END

!-------------------------------------------------------
! Important, first finished!!!!
!-------------------------------------------------------
SUBROUTINE SCGS
include 'incl.h'
integer IT,I,J
real*8 a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5

open(1,file="Residual.plt")
DO IT=1,MAXIT
! error residual set 0
      RESU=0.0
      RESV=0.0
      RESM=0.0
!block SOR
      DO J=1,NJ
      DO I=1,NI
!c for u(i  ,j)
!c BU means pressure gradient=dp*dx_i
		fe_u(1)=0.5d0*DY*(U(I,J)+U(I+1,J)) 
		fw_u(1)=0.5d0*DY*(U(I,J)+U(I-1,J))
		fn_u(1)=0.5d0*DX*(V(I,J)+V(I+1,J)) 
		fs_u(1)=0.5d0*DX*(V(I,J-1)+V(I+1,J-1))
		fe_u(0)=0.5d0*DY*(U(I-1,J)+U(I,J))
		fw_u(0)=0.5d0*DY*(U(I-1,J)+U(I-2,J))
		fn_u(0)=0.5d0*DX*(V(I-1,J)+V(I,J)) 
		fs_u(0)=0.5d0*DX*(V(I-1,J-1)+V(I,J-1))
		fe_v(1)=0.5d0*DY*(U(I,J)+U(I,J+1)) 
		fw_v(1)=0.5d0*DY*(U(I-1,J)+U(I-1,J+1))
		fn_v(1)=0.5d0*DX*(V(I,J+1)+V(I,J))
		fs_v(1)=0.5d0*DX*(V(I,J-1)+V(I,J))
		fe_v(0)=0.5d0*DY*(U(I,J-1)+U(I,J))
		fw_v(0)=0.5d0*DY*(U(I-1,J-1)+U(I-1,J))
		fn_v(0)=0.5d0*DX*(V(I,J)+V(I,J-1))
		fs_v(0)=0.5d0*DX*(V(I,J-2)+V(I,J-1))
if (Id_Escheme.eq.2)then
	call hybrid(I,J,b1,b2,b3,b4,b5)
else
	call solver(I,J,b1,b2,b3,b4,b5)
endif
!right above 2018.2.21
      a11=APU(0)
      a15=DY
      a22=APU(1)
      a25=-DY
      a33=APV(0)
      a35=DX
      a44=APV(1)
      a45=-DX
! under-relaxation
      a11=a11/URFU
      a22=a22/URFU
      a33=a33/URFV
      a44=a44/URFV

      a51=-DY
      a52= DY
      a53=-DX
      a54= DX
	  call BC_VANKA(I,J,a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5)

      CALL VANKA( &
     a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,&
     b1,b2,b3,b4,b5,&
     x1,x2,x3,x4,x5)
! correcting U,V and P
      U(I-1,J  )=U(I-1,J  )+x1
      U(I  ,J  )=U(I  ,J  )+x2
      V(I  ,J-1)=V(I  ,J-1)+x3
      V(I  ,J  )=V(I  ,J  )+x4
      P(I  ,J  )=P(I  ,J  )+x5
! calculate residuals at centroid of P-CV
      RESU=RESU+(abs(b1)+abs(b2))/2.
      RESV=RESV+(abs(b3)+abs(b4))/2.
      RESM=RESM+abs(b5)
      END DO
      END DO
! specify BCs for U and V
      CALL BCMOD
	  write(*,*)"iter=",IT,"error=",MAX(RESU,RESV,RESM)
	  write(1,*)IT,MAX(RESU,RESV,RESM)
      IF(MAX(RESU,RESV,RESM).LE.SORMAX.AND.IT.GT.1) GO TO 1000

END DO

 1000 CONTINUE
      RETURN
      END

!c-------------------------------------------------------
!c input:a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,
!c       &b1,b2,b3,b4,b5,
!c output:x1,x2,x3,x4,x5
!c-------------------------------------------------------
SUBROUTINE VANKA(&
     a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,&
     b1,b2,b3,b4,b5,&
     x1,x2,x3,x4,x5)
implicit none
real*8 DEN,r1,r2,r3,r4
real*8 a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5
      r1=a51/a11
      r2=a52/a22
      r3=a53/a33
      r4=a54/a44
      DEN=r1*a15+r2*a25+r3*a35+r4*a45
      x5=(r1*b1+r2*b2+r3*b3+r4*b4-b5)/DEN
      x1=(b1-a15*x5)/a11
      x2=(b2-a25*x5)/a22
      x3=(b3-a35*x5)/a33
      x4=(b4-a45*x5)/a44
RETURN
END

!c-------------------------------------------------------
!C! specify BCs for U and V
!driven velocity = 1m/s
!C******************************
!c-------------------------------------------------------
SUBROUTINE BCMOD
include 'incl.h'
	integer i,j
!BC up && down
	do i=0,NI
		U(i,NJ+1)=2.0-U(i,NJ)
		U(i,NJ+2)=2.0-U(i,NJ-1)
		U(i,0) =-U(i,1)
		U(i,-1)=-U(i,2)

		V(i,NJ)=0.0D0
		V(i,NJ+1)=-V(i,NJ-1)
		V(i,NJ+2)=-V(i,NJ-2)

		V(i,0)=0.0d0
		V(i,-1)=-V(i,1)
		V(i,-2)=-V(i,2)
	enddo
!left && right 
	do j=0,NJ
		U(0,j)=0.0d0
		U(-1,j)=-U(1,j)
		U(-2,j)=-U(2,j)

		U(NI,j)=0.0d0
		U(NI+1,j)=-U(NI-1,j)
		U(NI+2,j)=-U(NI-2,j)

		V(0,j)=-V(1,j)
		V(-1,j)=-V(2,j)

		V(NI+1,j)=-V(NI,j)
		V(NI+2,j)=-V(NI-1,j)
	enddo

      RETURN
      END


SUBROUTINE BC_VANKA(I,J,a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5)
include 'incl.h'
real*8		 a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5
integer I,J
!bc for continuity eq.
	if(I.eq.1) 	a51=0.0
	if(I.eq.NI)	a52=0.0
	if(J.eq.1)	a53=0.0
	if(J.eq.NJ)	a54=0.0
!bc for u
	if(I.eq.1) then
		a11=1.0
		a15=0.0
		b1=0.0
	endif
	if(I.eq.NI) then
		a22=1.0
		a25=0.0
		b2=0.0
	endif
!bc for v
	if(J.eq.1) then
		a33=1.0
		a35=0.0
		b3=0.0
	endif
	if(J.eq.NJ) then
		a44=1.0
		a45=0.0
		b4=0.0
	endif

return
end

subroutine output
      include 'incl.h'
	  integer I,J
	  real*8 XC,YC
	  real*8 UC(1:NI,1:NJ),VC(1:NI,1:NJ)
      DO I=1,NI
      DO J=1,NJ
      UC(I,J)=0.5*(U(I,J)+U(I-1,J))
      VC(I,J)=0.5*(V(I,J)+V(I,J-1))
      END DO
      END DO
!velocity field
      open(50,file='field.plt')
      WRITE(50,*)'VARIABLES = "x", "y", "u", "v", "p"'
      WRITE(50,*)'ZONE F=POINT, I=',NI, ', J=',NJ
      DO J=1,NJ
      DO I=1,NI
			XC=0.5*DX+(I-1)*DX
			YC=0.5*DY+(J-1)*DY
			WRITE(50,*)XC,YC,UC(I,J),VC(I,J),P(I,J)
      END DO
      END DO
      close(50)
!Center line for u 
		open(100,file="u_profile.plt")
		write(100,*)'variables="y","u"'
		do J=1,NJ
			YC=(J-0.5)*DY
			WRITE(100,*)YC,U(NI/2,J)
		enddo
		close(100)
!Center line for  v
		open(100,file="v_profile.plt")
		write(100,*)'variables="x","v"'
		do I=1,NI
			XC=(I-0.5)*DX
			WRITE(100,*)XC,V(I,NJ/2)
		enddo
		close(100)
return 
end













