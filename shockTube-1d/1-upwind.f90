!-------------------------------------------------------
	  subroutine upwind_1(Ka,Kb,fpa,fma,f)    !输入Kl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

        call flux_P_1upwind(Ka,Kb,fpa,fpx)
 	    call flux_M_1upwind(Ka,Kb,fma,fmx)
		f=fpx+fmx
	  return
	  end

!  1阶迎风格式(正通量)  周二晚上回来继续搞定
      subroutine flux_P_1upwind(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
	  real*8::w0,w1,q0,q1
	  w0=1.0d0
	  w1= 0.0d0
	  
	  q0=v(0)
	  q1=v(1)

	  v1=w0*q0+w1*q1

	 end
!-------------------------------------------------
!  1阶迎风格式(负通量)
      subroutine flux_M_1upwind(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb
	  real*8:: v(Ka:Kb),v1
	  real*8::w0,w1,w2,q0,q1,q2
	  w0=1.0d0
	  w1= 0.0d0
	  
	  q0=v(1)
	  q1=v(0)

	  v1=w0*q0+w1*q1
	
   end	  
    
!c---------------------------------------------------------------

