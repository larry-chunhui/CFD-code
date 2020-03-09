!-------------------------------------------------------
	  subroutine upwind_3(Ka,Kb,fpa,fma,f)    !输入Kl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

        call flux_P_3upwind(Ka,Kb,fpa,fpx)
 	    call flux_M_3upwind(Ka,Kb,fma,fmx)
		f=fpx+fmx
	  return
	  end
!-------------------------------------------------------
!  3阶迎风格式(正通量)
      subroutine flux_P_3upwind(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
	  real*8::w0,w1,w2,w3,q0,q1,q2,q3
	  w0=-1.0d0/6.0d0
	  w1= 5.0d0/6.0d0
	  w2= 2.0d0/6.0d0

	  
	  q0=v(-1)
	  q1=v(0)
	  q2=v(1)


	  v1=w0*q0+w1*q1+w2*q2

	 end
!-------------------------------------------------
!  3阶迎风格式(负通量)
      subroutine flux_M_3upwind(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
	  real*8::w0,w1,w2,w3,q0,q1,q2,q3
	  w0= 2.0d0/6.0d0
	  w1= 5.0d0/6.0d0
	  w2=-1.0d0/6.0d0

	  
	  q0=v(2)
	  q1=v(1)
	  q2=v(0)


	  v1=w0*q0+w1*q1+w2*q2


	
   end	  
    
!c---------------------------------------------------------------

