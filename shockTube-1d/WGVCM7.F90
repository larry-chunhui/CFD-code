! WGVC 格式（不含保单调限制）  ! 格式by 何志伟
!================================================================================
	  subroutine wgvc7(Ka,Kb,fpa,fma,f)    !输入Kl,Kr
!	  use var
	  implicit none
	  integer Ka,Kb
!	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

        call flux_P_wgvc7(Ka,Kb,fpa,fpx)
 	    call flux_M_wgvc7(Ka,Kb,fma,fmx)
		f=fpx+fmx
	  return
	  end
!		subroutine WGVC_M7_PO(v,mid_nf)
        subroutine flux_P_wgvc7(Ka,Kb,v,mid_nf)
			implicit none
			integer:: Ka,Kb
			real*8,dimension(Ka:Kb)::v
			real*8::mid_nf
			!---WENO-JS格式参数---
			real*8,dimension(0:3)::beta,d,w,alpha
			real*8,parameter::epsilon=1.d-40
			integer,parameter::p=2
			!---中间参数---
			integer::k
			real*8::h0,h1,h2,h3
			!---gvc----
			real*8::gm,gs,wgm,wgs,tao,sita
!           MP
   			real*8::minmod,minmod1
			real*8::d1,d2,d3,ul,md,lc,fmax,fmin
			real*8::kappa=2.d0



			d(0)=1.d0/35.d0
			d(1)=12.d0/35.d0
			d(2)=18.d0/35.d0
			d(3)=4.d0/35.d0

			h1=-2.d0/6.d0*v(-3)+9.d0/6.d0*v(-2)-18.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
			h2=-v(-3)+4.d0*v(-2)-5.d0*v(-1)+2.d0*v(0)
			h3=-v(-3)+3.d0*v(-2)-3.d0*v(-1)+v(0)
			beta(0)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3 

			h1=1.d0/6.d0*v(-2)-6.d0/6.d0*v(-1)+3.d0/6.d0*v(0)+2.d0/6.d0*v(1)
			h2=v(-1)-2.d0*v(0)+v(1) 
			h3=-v(-2)+3.d0*v(-1)-3.d0*v(0)+v(1) 
			beta(1)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3            

			h1=-2.d0/6.d0*v(-1)-3.d0/6.d0*v(0)+6.d0/6.d0*v(1)-1.d0/6.d0*v(2) 
			h2=v(-1)-2.d0*v(0)+v(1)   
			h3=-v(-1)+3.d0*v(0)-3.d0*v(1)+v(2)  
			beta(2)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3           
    
			h1=-11.d0/6.d0*v(0)+18.d0/6.d0*v(1)-9.d0/6.d0*v(2)+2.d0/6.d0*v(3)
			h2=2.d0*v(0)-5.d0*v(1)+4.d0*v(2)-1.d0*v(3)         
			h3=-v(0)+3.d0*v(1)-3.d0*v(2)+v(3)  
			beta(3)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3                 
  
			h0=0.d0
			do k=0,3
				alpha(k)=d(k)/((epsilon+beta(k))**p)
				h0=h0+alpha(k)
			enddo
	
			do k=0,3
				w(k)=alpha(k)/h0
			enddo
	
			!---4阶差分格式的通量---
			h0=-3.d0/12.d0*v(-3)+13.d0/12.d0*v(-2)-23.d0/12.d0*v(-1)+25.d0/12.d0*v(0)
			h1=1.d0/12.d0*v(-2)-5.d0/12.d0*v(-1)+13.d0/12.d0*v(0)+3.d0/12.d0*v(1)
			h2=-1.d0/12.d0*v(-1)+7.d0/12.d0*v(0)+7.d0/12.d0*v(1)-1.d0/12.d0*v(2)
			h3=3.d0/12.d0*v(0)+13.d0/12.d0*v(1)-5.d0/12.d0*v(2)+1.d0/12.d0*v(3)
			!=======================================================================================
			gm=1000.d0/3087.d0;gs=2087.d0/3087.d0
			tao=abs(beta(0)-beta(3))
			alpha(0)=gm*(1.d0+(tao/(beta(0)+epsilon))**3)
			alpha(1)=gs*(1.d0+(tao/(beta(3)+epsilon))**3)

			beta(0)=alpha(0)+alpha(1)
			wgm=alpha(0)/beta(0);wgs=alpha(1)/beta(0);

			sita=1.d0-wgm*wgs/(gm*gs)
			sita=sita**100*(101.d0-100.d0*sita)

			w(0)=(1-sita)*(0.0882d0*wgm)+sita*w(0)
			w(1)=(1-sita)*(0.2d0+0.441d0*wgm)+sita*w(1)
			w(2)=(1-sita)*(0.6d0-0.2646d0*wgm)+sita*w(2)
			w(3)=(1-sita)*(0.2d0-0.2646d0*wgm)+sita*w(3)
 !========================================================================================

			mid_nf=w(0)*h0+w(1)*h1+w(2)*h2+w(3)*h3
 
			return
		   endsubroutine

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!		subroutine WGVC_M7_NE(v,mid_nf)
	        subroutine flux_M_wgvc7(Ka,Kb,v,mid_nf)
			implicit none
			integer:: Ka,Kb
			real*8,dimension(Ka:Kb)::v
			real*8::mid_nf
			!---parameters for WENO-JS---
			real*8,dimension(0:3)::beta,d,w,alpha
			real*8,parameter::epsilon=1.d-40
			integer,parameter::p=2
			integer::k
			real*8::h0,h1,h2,h3
			!---gvc----
			real*8::gm,gs,wgm,wgs,tao,sita
!          MP
			real*8::minmod,minmod1
			real*8::d1,d2,d3,ul,md,lc,fmax,fmin
			real*8::kappa=2.d0    ! kappa=1.d0 (more robust)


			d(0)=1.d0/35.d0
			d(1)=12.d0/35.d0
			d(2)=18.d0/35.d0
			d(3)=4.d0/35.d0

			h1=2.d0/6.d0*v(4)-9.d0/6.d0*v(3)+18.d0/6.d0*v(2)-11.d0/6.d0*v(1)
			h2=-v(4)+4.d0*v(3)-5.d0*v(2)+2.d0*v(1)
			h3=v(4)-3.d0*v(3)+3.d0*v(2)-v(1)
			beta(0)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3 

			h1=-1.d0/6.d0*v(3)+6.d0/6.d0*v(2)-3.d0/6.d0*v(1)-2.d0/6.d0*v(0)
			h2=v(2)-2.d0*v(1)+v(0) 
			h3=v(3)-3.d0*v(2)+3.d0*v(1)-v(0) 
			beta(1)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3            

			h1=2.d0/6.d0*v(2)+3.d0/6.d0*v(1)-6.d0/6.d0*v(0)+1.d0/6.d0*v(-1) 
			h2=v(2)-2.d0*v(1)+v(0)
			h3=v(2)-3.d0*v(1)+3.d0*v(0)-v(-1)  
			beta(2)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3           
    
			h1=11.d0/6.d0*v(1)-18.d0/6.d0*v(0)+9.d0/6.d0*v(-1)-2.d0/6.d0*v(-2)
			h2=2.d0*v(1)-5.d0*v(0)+4.d0*v(-1)-1.d0*v(-2)         
			h3=v(1)-3.d0*v(0)+3.d0*v(-1)-v(-2)  
			beta(3)=h1*h1+13.d0/12.d0*h2*h2+1043.d0/960.d0*h3*h3+1.d0/12.d0*h1*h3            
  
			h0=0.d0
			do k=0,3
			alpha(k)=d(k)/((epsilon+beta(k))**p)
			h0=h0+alpha(k)
			enddo
	
			do k=0,3
			w(k)=alpha(k)/h0
			enddo
	
			!---flux of 4th order schemes--
			h0=-3.d0/12.d0*v(4)+13.d0/12.d0*v(3)-23.d0/12.d0*v(2)+25.d0/12.d0*v(1)
			h1=1.d0/12.d0*v(3)-5.d0/12.d0*v(2)+13.d0/12.d0*v(1)+3.d0/12.d0*v(0)
			h2=-1.d0/12.d0*v(2)+7.d0/12.d0*v(1)+7.d0/12.d0*v(0)-1.d0/12.d0*v(-1)
			h3=3.d0/12.d0*v(1)+13.d0/12.d0*v(0)-5.d0/12.d0*v(-1)+1.d0/12.d0*v(-2)
			!============================================
			gm=1000.d0/3087.d0;gs=2087.d0/3087.d0
			tao=abs(beta(0)-beta(3))
			alpha(0)=gm*(1.d0+(tao/(beta(0)+epsilon))**3)
			alpha(1)=gs*(1.d0+(tao/(beta(3)+epsilon))**3)

			beta(0)=alpha(0)+alpha(1)
			wgm=alpha(0)/beta(0);wgs=alpha(1)/beta(0);

			sita=1.d0-wgm*wgs/(gm*gs)
			sita=sita**100*(101.d0-100.d0*sita)

			w(0)=(1-sita)*(0.0882d0*wgm)+sita*w(0)
			w(1)=(1-sita)*(0.2d0+0.441d0*wgm)+sita*w(1)
			w(2)=(1-sita)*(0.6d0-0.2646d0*wgm)+sita*w(2)
			w(3)=(1-sita)*(0.2d0-0.2646d0*wgm)+sita*w(3)
			!============================================

			mid_nf=w(0)*h0+w(1)*h1+w(2)*h2+w(3)*h3

	
			return
		endsubroutine
	
! ==========================================

