! TENO5
!eps=1.e-10 no problem
!eps=1.e-40 the same as fulin-c++, problems
!-------------------------------------------------------

	  subroutine teno5_test(Ka,Kb,fpa,fma,f)    !输入Kl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

        call flux_P_teno5_test(Ka,Kb,fpa,fpx)
 	    call flux_M_teno5_test(Ka,Kb,fma,fmx)
		f=fpx+fmx
	  return
	  end
!!-------------------------------------------------------
!  5阶TENO格式（正通量）
      subroutine flux_P_teno5_test(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1,tau5
      real*8:: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23
	  real*8:: b1,b2,b0
      ep=1.d-10
	  C03=3.d0/10.d0
   	  C13=3.d0/5.d0
	  C23=1.d0/10.d0


         S0=13.d0/12.d0*(v(0)-2.d0*v(1)+v(2))**2+ 1.d0/4.d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
         S1=13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+ 1.d0/4.d0*(v(-1)-v(1))**2
         S2=13.d0/12.d0*(v(-2)-2.d0*v(-1)+v(0))**2+ 1.d0/4.d0*(v(-2)-4.d0*v(-1)+3.*v(0))**2

		 tau5=dabs(S2-S0)

         a0=(1+tau5/(S0+ep))**6
  	     a1=(1+tau5/(S1+ep))**6
	     a2=(1+tau5/(S2+ep))**6

		 b0=a0/(a0+a1+a2)
		 b1=a1/(a0+a1+a2)
		 b2=a2/(a0+a1+a2)
		 if(b0.ge.(1e-5)) then
			b0=1.0
		 else
		    b0=0.0
!			write(*,*)"a0=",a0
!			write(*,*)"a1=",a1
!			write(*,*)"a2=",a2
!			write(*,*)S0,S1,S2,tau5
!			stop
		 endif

		 if(b1.ge.(1e-5)) then
			b1=1.0
		 else
		    b1=0.0
		 endif

		 if(b2.ge.(1e-5)) then
			b2=1.0
		 else
		    b2=0.0
		 endif
!sequence different from fu code
		 a0=0.3*b0
		 a1=0.6*b1
		 a2=0.1*b2

 	     W0=a0/(a0+a1+a2)
         W1=a1/(a0+a1+a2)
	     W2=a2/(a0+a1+a2)
!write(*,*) a0,a1,a2
!write(*,*) W0,W1,W2

!	     q03=1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
!	     q13=-1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)
!	     q23=1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
!	     v1=W0*q03+W1*q13+W2*q23

	     q03=   1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)-v(0)
	     q13= -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)-v(0)
	     q23=1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)-v(0)
	     v1=v(0)+W0*q03+W1*q13+W2*q23
	 end
!-------------------------------------------------
!  5阶TENO格式 （负通量）
      subroutine flux_M_teno5_test(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
      real*8:: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23
	  real*8:: b1,b2,b0,tau5
	 ep=1.d-10
	 
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

        S0=13.d0/12.d0*(v(1)-2.d0*v(0)+v(-1))**2+  1.d0/4.d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
        S1=13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+  1.d0/4.d0*(v(2)-v(0))**2
        S2=13.d0/12.d0*(v(3)-2.d0*v(2)+v(1))**2+ 1.d0/4.d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
!Chunhui change
		 tau5=dabs(S2-S0)

         a0=(1+tau5/(S0+ep))**6
  	     a1=(1+tau5/(S1+ep))**6
	     a2=(1+tau5/(S2+ep))**6
		 b0=a0/(a0+a1+a2)
		 b1=a1/(a0+a1+a2)
		 b2=a2/(a0+a1+a2)
		 if(b0.ge.(1e-5)) then
			b0=1.0
		 else
		    b0=0.0
		 endif

		 if(b1.ge.(1e-5)) then
			b1=1.0
		 else
		    b1=0.0
		 endif

		 if(b2.ge.(1e-5)) then
			b2=1.0
		 else
		    b2=0.0
		 endif
!sequence different from fu code
		 a0=0.3*b0
		 a1=0.6*b1
		 a2=0.1*b2


 	    W0=a0/(a0+a1+a2)
        W1=a1/(a0+a1+a2)
	    W2=a2/(a0+a1+a2)
	    q03=1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)-v(1)
	    q13=-1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)-v(1)
	    q23=1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)-v(1)
 	
	    v1=v(1)+W0*q03+W1*q13+W2*q23
	
   end	  
    
!c---------------------------------------------------------------

