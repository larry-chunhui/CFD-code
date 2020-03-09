subroutine center6(Ka,Kb,fpa,fma,f,delta_p,delta_m)    ! ‰»ÎKl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb,k
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8::f,f_p,f_m,d1,d2,d3,delta_p,delta_m
	  real*8 fpa(Ka:Kb),fma(Ka:Kb),v(Ka:Kb),v_p(Ka:Kb),v_m(Ka:Kb)
      v_p=fpa
      v_m=fma
      d1=3.0d0/4.0d0
      d2=-3.0d0/20.0d0
      d3=1.0d0/60.0d0
      k=0
      f_p=delta_p*((d1+d2+d3)*(v_p(k)+v_p(k+1))+(d2+d3)*(v_p(k-1)+v_p(k+2))+d3*(v_p(k-2)+v_p(k+3)))
      f_m=delta_m*((d1+d2+d3)*(v_m(k)+v_m(k+1))+(d2+d3)*(v_m(k-1)+v_m(k+2))+d3*(v_m(k-2)+v_m(k+3))) 
      f=f_p+f_m
	  return
    end
    