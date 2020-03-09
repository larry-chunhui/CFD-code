!-------------------------------------------------------

	  subroutine weno5(Ka,Kb,fpa,fma,f)    !����Kl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

        call flux_P_weno5(Ka,Kb,fpa,fpx)
 	    call flux_M_weno5(Ka,Kb,fma,fmx)
		f=fpx+fmx
	  return
	  end
!!-------------------------------------------------------
!  5��WENO��ʽ����ͨ����
      subroutine flux_P_weno5(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
      real*8:: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23
      ep=1.d-6
	  C03=3.d0/10.d0
   	  C13=3.d0/5.d0
	  C23=1.d0/10.d0


         S0=13.d0/12.d0*(v(0)-2.d0*v(1)+v(2))**2+ 1.d0/4.d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
         S1=13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+ 1.d0/4.d0*(v(-1)-v(1))**2
         S2=13.d0/12.d0*(v(-2)-2.d0*v(-1)+v(0))**2+ 1.d0/4.d0*(v(-2)-4.d0*v(-1)+3.*v(0))**2
         a0=C03/((ep+S0)**2)
  	     a1=C13/((ep+S1)**2)
	     a2=C23/((ep+S2)**2)
 	     W0=a0/(a0+a1+a2)
         W1=a1/(a0+a1+a2)
	     W2=a2/(a0+a1+a2)
	     q03=1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
	     q13=-1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)
	     q23=1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
	     v1=W0*q03+W1*q13+W2*q23

	 end
!-------------------------------------------------
!  5��WENO��ʽ ����ͨ����
      subroutine flux_M_weno5(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
      real*8:: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23

	 ep=1.d-6
	 
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

        S0=13.d0/12.d0*(v(1)-2.d0*v(0)+v(-1))**2+  1.d0/4.d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
        S1=13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+  1.d0/4.d0*(v(2)-v(0))**2
        S2=13.d0/12.d0*(v(3)-2.d0*v(2)+v(1))**2+ 1.d0/4.d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
        a0=C03/((ep+S0)**2)
	    a1=C13/((ep+S1)**2)
	    a2=C23/((ep+S2)**2)
 	    W0=a0/(a0+a1+a2)
        W1=a1/(a0+a1+a2)
	    W2=a2/(a0+a1+a2)
	    q03=1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)
	    q13=-1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)
	    q23=1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)
 	
	    v1=W0*q03+W1*q13+W2*q23
	
   end	  
    
!c---------------------------------------------------------------
