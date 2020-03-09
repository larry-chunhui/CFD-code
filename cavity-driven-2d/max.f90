program main
implicit none
real*8 a,b,c,d,e
! c hello everyone
	a=0.0d0
	b=1.0d0
	c=-7.0d0
	d=9.0d0
	e=max(a,b)
	write(*,*)"a=",a,"b=",b,"max=",e
	e=max(a,c)
	write(*,*)"a=",a,"c=",c,"max=",e
	e=max(b,c)
	write(*,*)"b=",b,"c=",c,"max=",e
	e=a+b+ 	&
          c+d
	write(*,*)"-----------------------------"
	write(*,*)"e=3?",e
stop
end
