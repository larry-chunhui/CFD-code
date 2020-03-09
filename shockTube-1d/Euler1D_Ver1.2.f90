!--------------------------------------------------------------------------------------------
! Euler 1D, Platform of testing schemes
!wen5,wen7,wen5z,1upwind,3upwind, wgvc7scheme all support. 2018.04.02
!---------------------------------------------------------------------------------
   module var
   implicit none
   integer:: nx,Istep                                  ! �������,ʱ�䲽
   integer:: Iflag_splitting, Iflag_Characteric, Iflag_init, Iflag_scheme !,Iflag_smooth
   integer:: KL,KR  ! ����j+1/2ʹ�õ�������ܵ���±�  (���磬 ʹ��j-3,j-2,...j+4 �����j+1/2����ֵ�� �� KL=-3, KR=4)
   integer,parameter:: Split_Steger_Warming=0, Split_LF_local=1, Split_LF_Global=2, Split_VanLeer=3
   integer,parameter:: Init_Sod=0, Init_Shu_Osher=1, Init_lien=2,Iflag_smooth=1
   real*8,parameter:: gamma=1.4d0
   real*8,dimension(:),Pointer:: d,u,p,c          ! ԭʼ����
   real*8,dimension(:,:),Pointer:: fp,fm,flux         ! ������ͨ�� (4,nx) , ͨ��
   real*8,dimension(:,:),Pointer:: Q,Qn           ! �غ���� Q(4,nx) , Qn ��һʱ�䲽��Q
   real*8,dimension(:),Pointer:: x    ! ����
   real*8:: Ma, Cp, Cv, dt, dx, tt,t_end
   end module var

!------------------------------------------------------------------------------
   program Euler1D
   use var
   implicit none
   integer:: KRK,i
   real t1,t2
   
   call cpu_time(t1)
   call init              ! ��ʼ��
  !                        ʱ���ƽ�  
write(*,*)"1e-5",1e-5  
    do    ! ʱ���ƽ�
     Qn=Q
     do KRK=1,3
	   call comput_original_variables        ! ����ԭʼ���� (d,u,p)
       call splitting                        ! ͨ������
 !            ����ͨ��
		if( Iflag_Characteric .eq. 0) then
	      call comput_flux                    ! ����ͨ������ʹ�þֲ������ֽ⣩                      
        else 
	      call comput_flux_Characteric        ! ����ͨ�� ��ʹ�þֲ������ֽ⣩                
	    endif

	   call time_advance (KRK) 
     enddo
     print*, "Istep, tt,u=", Istep,tt,u(nx/2)
     if(tt .gt. t_end) exit
    enddo
    call cpu_time(t2)
    open(1,file="cputime.dat")
	write(1,*)"cpu time=",t2-t1
    close(1)

    call comput_original_variables        ! ����ԭʼ���� (d,u,p)
    open(99,file="sw-1d.plt")
	write(99,*)"variables=x,d,u,M,p"
	do i=1,nx
	 write(99,"(4f16.8)") x(i),d(i),u(i),u(i)/c(i),p(i)
	enddo
    close(99)
   end
!---------------------------------------------------------------------------------
   subroutine init
   use var
   implicit none
   integer:: i
   real*8:: d1,u1,p1,d2,u2,p2
!----------------------------------------------    
    open(99,file="control.in")
	read(99,*)
	read(99,*) nx, Iflag_splitting,	Iflag_Characteric, Iflag_init,Iflag_scheme
	read(99,*)
	read(99,*) 	KL,KR, dt,t_end
!---------------------------------------------
    Istep=0
    tt=0.d0

    allocate(x(nx),d(nx),u(nx),p(nx),c(nx),fp(4,nx),fm(4,nx),flux(4,nx),Q(4,nx),Qn(4,nx))
 
    if(Iflag_init .eq. Init_Sod) then   ! Sod probrem
   	  dx=1.d0/(nx-1.d0)
      do i=1,nx
	  x(i)=(i-1.d0)*dx
	  enddo
	  d1=1.d0;  u1=0.d0 ;    p1=1.d0
      d2=0.1d0; u2=0.d0 ;    p2=0.125d0
      do  i=1,nx        
      if(i.lt.nx/2) then
      Q(1,i)=d1
      Q(2,i)=d1*u1
      Q(3,i)=p1/(gamma-1.d0)+d1*u1*u1*0.5d0
	  else
      Q(1,i)=d2
      Q(2,i)=d2*u2
      Q(3,i)=p2/(gamma-1.d0)+d2*u2*u2*0.5d0
      endif
	 enddo
    
	else if(Iflag_init .eq. Init_Shu_Osher) then  ! Shu_Osher problem
   	   dx=10.d0/(nx-1.d0)
        do i=1,nx
	     x(i)=(i-1.d0)*dx
 	    enddo
	   do  i=1,nx        
	   if(x(i).le. 1.d0) then
         d(i)=3.857143d0
	     p(i)=10.33333d0
	     u(i)=2.629369d0
	   else
         d(i)=1.d0+0.2*sin(5.d0*x(i))
	     p(i)=1.d0
	     u(i)=0.d0
       endif
       Q(1,i)=d(i)
       Q(2,i)=d(i)*u(i)
       Q(3,i)=p(i)/(gamma-1.d0)+d(i)*u(i)*u(i)*0.5d0
      enddo

     else if(Iflag_init .eq. Init_lien) then   ! Sod probrem
   	  dx=20.0d0/(nx-1.d0)
      do i=1,nx
	  x(i)=(i-1.d0)*dx-10.0d0
	  enddo
	  d1=1.d0;  u1=0.d0 ;    p1=100000.d0
      d2=0.125d0; u2=0.d0 ;    p2=10000.d0
      do  i=1,nx        
      if(i.lt.nx/2) then
      Q(1,i)=d1
      Q(2,i)=d1*u1
      Q(3,i)=p1/(gamma-1.d0)+d1*u1*u1*0.5d0
	  else
      Q(1,i)=d2
      Q(2,i)=d2*u2
      Q(3,i)=p2/(gamma-1.d0)+d2*u2*u2*0.5d0
      endif
	 enddo
    endif
   end subroutine init



    subroutine time_advance(KRK)
    use var
    implicit none
    integer:: KRK,i,m
    real*8:: du
	 
	  do i=2-KL,nx-KR
      do m=1,3 
	    du=-(flux(m,i)-flux(m,i-1))/dx
	    
	   if(KRK.eq.1) then
        Q(m,i)=Qn(m,i)+dt*du
       else if (KRK .eq. 2) then
	    Q(m,i)=3.d0/4.d0*Qn(m,i)+1.d0/4.d0*(Q(m,i)+dt*du)
	   else if (KRK .eq. 3) then
	    Q(m,i)=1.d0/3.d0*Qn(m,i)+2.d0/3.d0*(Q(m,i)+dt*du)
       endif

	  enddo
      enddo
   
     if(KRK .eq. 3) then
	   tt=tt+dt
	   Istep=Istep+1
	 endif
   end
 

!----------------------------------------------------------------------------------
! �����غ���������ԭʼ����(d,u,p)
   subroutine comput_original_variables
    use var
	implicit none
    integer:: i
	do i=1,nx
      d(i)=Q(1,i)                                       ! �ܶ�
      u(i)=Q(2,i)/d(i)                                  ! �ٶ�
      p(i)=(gamma-1.d0)*(Q(3,i)-0.5d0*Q(2,i)*u(i))      ! ѹ��
      c(i)=sqrt(gamma*p(i)/d(i))                        ! ����
	enddo
   end subroutine comput_original_variables

! -----ͨ������----------------------------------------------------------------------
   subroutine splitting
    use var
	implicit none
	if(Iflag_splitting .eq. Split_LF_local) then                ! ����L-F���� (�����������ֵ����)
	   call Splitting_LF_local
    else if (Iflag_splitting .eq. Split_LF_Global) then         ! ȫ��L-F���� (ȫ���������ֵ����)
	   call Splitting_LF_Global  
	
	else if(Iflag_splitting .eq. Split_Steger_Warming) then
      call Splitting_Steger_Warming
	else if (Iflag_splitting .eq. Split_VanLeer) then
	  call Splitting_VanLeer
	endif
	end
!-----------------------------------------------------------------------------------------
!------����L-F ����
   subroutine Splitting_LF_local
    use var
	implicit none
    integer:: i
	real*8:: lmax
	do i=1,nx
     lmax=dmax1(abs(u(i)),abs(u(i)-c(i)),abs(u(i)+c(i)))        ! �����������ֵ
     fp(1,i)=0.5d0*(d(i)*u(i)+lmax*Q(1,i))
     fp(2,i)=0.5d0*(d(i)*u(i)*u(i)+p(i)+lmax*Q(2,i))
	 fp(3,i)=0.5d0*((Q(3,i)+p(i))*u(i)+lmax*Q(3,i))
     fm(1,i)=0.5d0*(d(i)*u(i)-lmax*Q(1,i))
     fm(2,i)=0.5d0*(d(i)*u(i)*u(i)+p(i)-lmax*Q(2,i))
	 fm(3,i)=0.5d0*((Q(3,i)+p(i))*u(i)-lmax*Q(3,i))
    enddo
   end subroutine Splitting_LF_local
!-------ȫ��L-F����   
   subroutine Splitting_LF_global
    use var
	implicit none
    integer:: i
	real*8:: lmax,lmax1
	 lmax=0.d0 
 	 do i=1,nx
	   lmax1=dmax1(abs(u(i)),abs(u(i)+c(i)),abs(u(i)-c(i)))        ! �����������ֵ
       lmax=dmax1(lmax1,lmax)
	 enddo

	do i=1,nx
     fp(1,i)=0.5d0*(d(i)*u(i)+lmax*Q(1,i))
     fp(2,i)=0.5d0*(d(i)*u(i)*u(i)+p(i)+lmax*Q(2,i))
	 fp(3,i)=0.5d0*((Q(3,i)+p(i))*u(i)+lmax*Q(3,i))
     fm(1,i)=0.5d0*(d(i)*u(i)-lmax*Q(1,i))
     fm(2,i)=0.5d0*(d(i)*u(i)*u(i)+p(i)-lmax*Q(2,i))
	 fm(3,i)=0.5d0*((Q(3,i)+p(i))*u(i)-lmax*Q(3,i))
    enddo
   end subroutine Splitting_LF_global
!---------------------------------------------------------
  	  
    subroutine Splitting_steger_warming
      use var
      implicit none
      real*8:: El1,El2,El3,El1p,El1m,El2p,El2m,El3p,El3m,tmp,WP,WM
      integer:: i

      do i=1,nx 
	    El1=u(i)
        El2=u(i)-c(i)
        El3=u(i)+c(i)
        El1p=0.5d0*(El1+abs(El1))
        El1M=0.5d0*(El1-abs(El1))
        El2p=0.5d0*(El2+abs(El2))
        El2M=0.5d0*(El2-abs(El2))
        El3p=0.5d0*(El3+abs(El3))
        El3M=0.5d0*(El3-abs(El3))
        tmp=d(i)/(2.d0*gamma)

        fp(1,i)=tmp*(2.d0*(gamma-1.d0)*El1p+El2p+El3p) 
        fm(1,i)=tmp*(2.d0*(gamma-1.d0)*El1M+El2M+El3M) 
        fp(2,i)=tmp*(2.d0*(gamma-1.d0)*El1p*u(i)+El2p*(u(i)-c(i))+El3p*(u(i)+c(i)) )
        fm(2,i)=tmp*(2.d0*(gamma-1.d0)*El1M*u(i)+El2M*(u(i)-c(i))+El3M*(u(i)+c(i)) )
        WP=(3.d0-gamma)*(El2p+El3p)*c(i)*c(i)/(2.d0*(gamma-1.d0))
        WM=(3.d0-gamma)*(El2M+El3M)*c(i)*c(i)/(2.d0*(gamma-1.d0))
        fp(3,i)=tmp*((gamma-1.d0)*El1p*u(i)*u(i)+El2p/2.d0*(u(i)-c(i))**2 +Wp+El3p/2.d0*(u(i)+c(i))**2)
        fm(3,i)=tmp*((gamma-1.d0)*El1M*u(i)*u(i)+El2M/2.d0*(u(i)-c(i))**2 +WM+El3M/2.d0*(u(i)+c(i))**2)
     enddo

    end subroutine Splitting_steger_warming






!--------------------------------------------------------------
    subroutine Splitting_VanLeer
    use var
    implicit none
	 real*8:: Max,Mp,Mm
     integer:: i
	  do i=1,nx
        Max=u(i)/c(i)

         if(Max>=1.d0) then
          fp(1,i)=d(i)*u(i)
          fp(2,i)=d(i)*u(i)*u(i)+p(i)
          fp(3,i)=(Q(3,i)+p(i))*u(i)
          fm(1,i)=0.d0; fm(2,i)=0.d0; fm(3,i)=0.d0
         else if(Max <=-1.d0) then
          fm(1,i)=d(i)*u(i)
          fm(2,i)=d(i)*u(i)*u(i)+p(i)
          fm(3,i)=(Q(3,i)+p(i))*u(i)
          fp(1,i)=0.d0; fp(2,i)=0.d0; fp(3,i)=0.d0
		 else   
          Mp= 0.25d0*(Max+1.d0)**2*d(i)*c(i)
          Mm=-0.25d0*(Max-1.d0)**2*d(i)*c(i)
          fp(1,i)=Mp
		  fp(2,i)=Mp*((gamma-1.d0)*u(i)+2.d0*c(i))/gamma
		  fp(3,i)=Mp*((gamma-1.d0)*u(i)+2.d0*c(i))**2/(2.d0*(gamma*gamma-1.d0))
          fm(1,i)=Mm
		  fm(2,i)=Mm*((gamma-1.d0)*u(i)-2.d0*c(i))/gamma
 		  fm(3,i)=Mm*((gamma-1.d0)*u(i)-2.d0*c(i))**2/(2.d0*(gamma*gamma-1.d0))
        endif
      enddo

	end subroutine Splitting_VanLeer
!-------------------------------------------------------------

! ����ͨ��, ֱ�Ӽ��㣬��ʹ�þֲ������ֽ�

     subroutine comput_flux                      
      use var
      implicit none
	  integer:: m,j
	  real*8,dimension(:),allocatable:: fpa,fma
	  real*8:: fpx,fmx
	  allocate(fpa(KL:KR),fma(KL:KR))
      do m=1,3
	   do j=1-KL,nx-KR
	    fpa(:)=fp(m,j+KL:j+KR)
        fma(:)=fm(m,j+KL:j+KR)
        call flux_P_weno5_z(KL,KR,fpa,fpx)
 	    call flux_M_weno5_z(KL,KR,fma,fmx)

!		call scheme
        flux(m,j)=fpx+fmx
	   enddo
	 enddo
     deallocate(fpa,fma)
	end

!--------------------------------------------------------------------------------------------
    subroutine comput_flux_Characteric                      ! ����ͨ�� (���þֲ������ֽ�)
      use var
      implicit none
	  integer:: m,j,k,k1
	  real*8,dimension(:,:),allocatable:: fpa,fma          ! ��������ͨ��               
	  real*8:: d1,u1,H1,d2,u2,H2,r1,r2,r0,u0,H0,c0,t1,s11,s12,s13,s21,s22,s23,s31,s32,s33
	  real*8:: t2,D11,D12,D13,D21,D22,D23,D31,D32,D33
      real*8:: fa0(3),fx(3)

 !-------------------------------------------------------------     
	   allocate(fpa(KL:KR,3),fma(KL:KR,3))


! Roe average 
         do j=1-KL,nx-KR
!        ������j+1/2���������� ��ʹ��Roeƽ����  ! Ҳ��ʹ�ü򵥵�����ƽ��
        
	      d1=d(j)
          u1=u(j)
          H1=gamma/(gamma-1.d0)*p(j)/d(j)+u(j)*u(j)*0.5d0            ! ����
	      d2=d(j+1)
          u2=u(j+1)
          H2=gamma/(gamma-1.d0)*p(j+1)/d(j+1)+u(j+1)*u(j+1)*0.5d0

!        Roe ƽ��
          r1=sqrt(d1) ; r2= sqrt(d2); r0= r1+r1
          u0=(r1*u1+r2*u2)/r0
	      H0=(r1*H1+r2*H2)/r0                     
	      c0=sqrt((gamma-1.d0)*(H0-u0*u0*0.5d0))                  ! j+1/2�㴦������

!  ��������S  (A������������)
	      t1=(gamma-1.d0)/c0
          s11=u0*u0/2.d0-c0*c0/(gamma-1.d0)
	      s12=-u0
	      s13=1.d0
	      s21=-u0-t1*u0*u0/2.d0
	      s22=1.d0+t1*u0
	      s23=-t1
	      s31=-u0+t1*u0*u0/2.d0
	      s32=1.d0-t1*u0
	      s33=t1

! ��������S^-1  (S�������A������������)
	      t2=(gamma-1.d0)/(c0*c0)
		  D11=-t2
	      D12=-0.5d0/c0
	      D13=0.5d0/c0
	      D21=-t2*u0
	      D22=-(u0-c0)/(2.d0*c0)
	      D23=(u0+c0)/(2.d0*c0)
	      D31=-t2*u0*u0/2.d0
	      D32=-(u0*u0/2.d0+1.d0/t2-u0*c0)/(2.d0*c0)
	      D33=(u0*u0/2.d0+1.d0/t2+u0*c0)/(2.d0*c0)


 !  �ã��ң���������S����������ܵ��ϵ�ÿһ���������ֽ���ͨ�������任�������ռ�
		 do k=KL,KR                 ! ����f(j+1/2)ʹ�õ�������ܵ� (����j-2,j-1, ... j+3 ʱ�� KL=-2, KR=3) 
           k1=j+k		
!          ��ͨ�� ���任�������ռ䣩
      	   fpa(k,1)=s11*fp(1,k1)+s12*fp(2,k1)+s13*fp(3,k1)
	       fpa(k,2)=s21*fp(1,k1)+s22*fp(2,k1)+s23*fp(3,k1)
	       fpa(k,3)=s31*fp(1,k1)+s32*fp(2,k1)+s33*fp(3,k1)
!          ��ͨ�� ���任�������ռ�)       	   
		   fma(k,1)=s11*fm(1,k1)+s12*fm(2,k1)+s13*fm(3,k1)
	       fma(k,2)=s21*fm(1,k1)+s22*fm(2,k1)+s23*fm(3,k1)
	       fma(k,3)=s31*fm(1,k1)+s32*fm(2,k1)+s33*fm(3,k1)
        
		 enddo

!    ���ò�ָ�ʽ������j+1/2����ͨ�� (��������ͨ��)  
         do m=1,3
!          call flux_P_weno5_z(KL,KR,fpa(:,m),fpx(m))
! 	      call flux_M_weno5_z(KL,KR,fma(:,m),fmx(m))
!		  call scheme(Ka,Kb,fpa,fma,f,Iflag_scheme)
		  call scheme(KL,KR,fpa(:,m),fma(:,m),fx(m),Iflag_scheme)
  	      fa0(m)=fx(m)
		 enddo

!    ʹ����任(D=S^-1), �õ�����ռ��ͨ�� ��j+1/2�㣩
          flux(1,j)=D11*fa0(1)+D12*fa0(2)+D13*fa0(3)
	      flux(2,j)=D21*fa0(1)+D22*fa0(2)+D23*fa0(3)
	      flux(3,j)=D31*fa0(1)+D32*fa0(2)+D33*fa0(3)
       enddo
     deallocate(fpa,fma)
	end
