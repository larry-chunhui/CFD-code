!-------------------------------------------------------

	  subroutine weno7(Ka,Kb,fpa,fma,f)    !输入Kl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

        call flux_P_weno7(Ka,Kb,fpa,fpx)
 	    call flux_M_weno7(Ka,Kb,fma,fmx)
		f=fpx+fmx
	  return
	  end
!____________________________________________________
! 7th order WENO scheme ( Jiang & Shu )
      subroutine flux_P_weno7(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb,k
	  real*8:: v(Ka:Kb),v1
      real*8:: ep,C0,C1,C2,C3,S10,S11,S12,S13,S20,S21,S22,S23,S30,S31,S32,S33,S0,S1,S2,S3
	  real*8:: a0,a1,a2,a3,W0,W1,W2,W3,q0,q1,q2,q3

       ep=1.d-7
	   C0=1.d0/35.d0
   	   C1=12.d0/35.d0
	   C2=18.d0/35.d0
	   C3=4.d0/35.d0

       k=0
         S10=-2.d0/6.d0*v(k-3)+9.d0/6.d0*v(k-2)-18.d0/6.d0*v(k-1)  +11.d0/6.d0*v(k)   ! 1 阶导数
         S11=1.d0/6.d0*v(k-2)-6.d0/6.d0*v(k-1)+3.d0/6.d0*v(k)      +2.d0/6.d0*v(k+1)
         S12=-2.d0/6.d0*v(k-1)-3.d0/6.d0*v(k)+6.d0/6.d0*v(k+1)     -1.d0/6.d0*v(k+2)
         S13=-11.d0/6.d0*v(k)+18.d0/6.d0*v(k+1)-9.d0/6.d0*v(k+2)   +2.d0/6.d0*v(k+3)

         S20=-v(k-3)+4.d0*v(k-2)-5.d0*v(k-1)+2.d0*v(k)              ! 2 阶导数
         S21=             v(k-1)-2.d0*v(k)  +v(k+1)         
         S22=             v(k)  -2.d0*v(k+1)+v(k+2)         
         S23=2.d0*v(k)-5.d0*v(k+1)+4.d0*v(k+2)-1.d0*v(k+3)         

         S30=-v(k-3)+3.d0*v(k-2)-3.d0*v(k-1)+v(k)                   ! 3 阶导数                
         S31=-v(k-2)+3.d0*v(k-1)-3.d0*v(k  )+v(k+1)                 
         S32=-v(k-1)+3.d0*v(k  )-3.d0*v(k+1)+v(k+2)                 
         S33=-v(k  )+3.d0*v(k+1)-3.d0*v(k+2)+v(k+3)                 

         S0=S10*S10+13.d0/12.d0*S20*S20 +1043.d0/960.d0*S30*S30 +1.d0/12.d0*S10*S30
         S1=S11*S11+13.d0/12.d0*S21*S21 +1043.d0/960.d0*S31*S31 +1.d0/12.d0*S11*S31
         S2=S12*S12+13.d0/12.d0*S22*S22 +1043.d0/960.d0*S32*S32 +1.d0/12.d0*S12*S32
         S3=S13*S13+13.d0/12.d0*S23*S23 +1043.d0/960.d0*S33*S33 +1.d0/12.d0*S13*S33

        a0=C0/((ep+S0)**2) ;  a1=C1/((ep+S1)**2) ;   a2=C2/((ep+S2)**2) ; a3=C3/((ep+S3)**2)
 	    W0=a0/(a0+a1+a2+a3) ;    W1=a1/(a0+a1+a2+a3) ;   W2=a2/(a0+a1+a2+a3) ;    W3=a3/(a0+a1+a2+a3)

!  4阶差分格式的通量
	   q0=-3.d0/12.d0*v(k-3)+13.d0/12.d0*v(k-2)-23.d0/12.d0*v(k-1)   +25.d0/12.d0*v(k)
	   q1=1.d0/12.d0*v(k-2)-5.d0/12.d0*v(k-1)+13.d0/12.d0*v(k)   +3.d0/12.d0*v(k+1)
	   q2=-1.d0/12.d0*v(k-1)+7.d0/12.d0*v(k)+7.d0/12.d0*v(k+1)  -1.d0/12.d0*v(k+2)
	   q3=3.d0/12.d0*v(k)+13.d0/12.d0*v(k+1)-5.d0/12.d0*v(k+2)  +1.d0/12.d0*v(k+3)

!  由4个4阶差分格式组合成1个7阶差分格式
	   v1=W0*q0+W1*q1+W2*q2+W3*q3
	 
   	  end subroutine flux_P_weno7

!C-------------------------------------------------
!C----- 对于负通量：
      subroutine flux_M_weno7(Ka,Kb,v,v1)
      implicit none
	  integer:: Ka,Kb,k
	  real*8:: v(Ka:Kb),v1
      real*8:: ep,C0,C1,C2,C3,S10,S11,S12,S13,S20,S21,S22,S23,S30,S31,S32,S33,S0,S1,S2,S3
	  real*8:: a0,a1,a2,a3,W0,W1,W2,W3,q0,q1,q2,q3

       ep=1.d-7
	   C0=1.d0/35.d0
   	   C1=12.d0/35.d0
	   C2=18.d0/35.d0
	   C3=4.d0/35.d0

      k=1
         S10=-2.d0/6.d0*v(k+3)+9.d0/6.d0*v(k+2)-18.d0/6.d0*v(k+1)   +11.d0/6.d0*v(k)   ! 1 阶导数
         S11=1.d0/6.d0*v(k+2)-6.d0/6.d0*v(k+1)+3.d0/6.d0*v(k)       +2.d0/6.d0*v(k-1)
         S12=-2.d0/6.d0*v(k+1)-3.d0/6.d0*v(k)+6.d0/6.d0*v(k-1)      -1.d0/6.d0*v(k-2)
         S13=-11.d0/6.d0*v(k)+18.d0/6.d0*v(k-1)-9.d0/6.d0*v(k-2)    +2.d0/6.d0*v(k-3)

         S20=-v(k+3)+4.d0*v(k+2)-5.d0*v(k+1)+2.d0*v(k)              ! 2 阶导数
         S21=             v(k+1)-2.d0*v(k)  +v(k-1)         
         S22=             v(k)  -2.d0*v(k-1)+v(k-2)         
         S23=2.d0*v(k)-5.d0*v(k-1)+4.d0*v(k-2)-1.d0*v(k-3)         

         S30=-v(k+3)+3.d0*v(k+2)-3.d0*v(k+1)+v(k)                   ! 3 阶导数                
         S31=-v(k+2)+3.d0*v(k+1)-3.d0*v(k  )+v(k-1)                 
         S32=-v(k+1)+3.d0*v(k  )-3.d0*v(k-1)+v(k-2)                 
         S33=-v(k  )+3.d0*v(k-1)-3.d0*v(k-2)+v(k-3)                 

         S0=S10*S10+13.d0/12.d0*S20*S20  +1043.d0/960.d0*S30*S30 +1.d0/12.d0*S10*S30
         S1=S11*S11+13.d0/12.d0*S21*S21  +1043.d0/960.d0*S31*S31 +1.d0/12.d0*S11*S31
         S2=S12*S12+13.d0/12.d0*S22*S22  +1043.d0/960.d0*S32*S32 +1.d0/12.d0*S12*S32
         S3=S13*S13+13.d0/12.d0*S23*S23  +1043.d0/960.d0*S33*S33 +1.d0/12.d0*S13*S33

       a0=C0/((ep+S0)**2) ;  a1=C1/((ep+S1)**2) ;   a2=C2/((ep+S2)**2) ;   a3=C3/((ep+S3)**2)
	   W0=a0/(a0+a1+a2+a3);  W1=a1/(a0+a1+a2+a3);  W2=a2/(a0+a1+a2+a3) ;  W3=a3/(a0+a1+a2+a3)

!  4阶差分格式的通量
	 q0=-3.d0/12.d0*v(k+3)+13.d0/12.d0*v(k+2)-23.d0/12.d0*v(k+1)  +25.d0/12.d0*v(k)
	 q1=1.d0/12.d0*v(k+2)-5.d0/12.d0*v(k+1)+13.d0/12.d0*v(k)      +3.d0/12.d0*v(k-1)
	 q2=-1.d0/12.d0*v(k+1)+7.d0/12.d0*v(k)+7.d0/12.d0*v(k-1)     -1.d0/12.d0*v(k-2)
	 q3=3.d0/12.d0*v(k)+13.d0/12.d0*v(k-1)-5.d0/12.d0*v(k-2)     +1.d0/12.d0*v(k-3)

!  由4个4阶差分格式组合成1个7阶差分格式
	   v1=W0*q0+W1*q1+W2*q2+W3*q3
	 
	 end  subroutine flux_M_weno7


















