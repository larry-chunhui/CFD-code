subroutine scheme(Ka,Kb,fpa,fma,f,Iflag_scheme)       !明天上午再来搞定,周一下午再来搞定
	  implicit none
	  integer Ka,Kb,Iflag_scheme

	  real*8:: fpx,fmx,f
	  real*8 fpa(Ka:Kb),fma(Ka:Kb)

if(Iflag_scheme.eq.0) then
	call weno5_z(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.1) then
	call weno5(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.2) then
	call weno7(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.3) then
	call upwind_1(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.4) then
	call upwind_3(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.5) then
	call wgvc7(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.6) then
	call teno5(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.61) then
	call teno5_test(Ka,Kb,fpa,fma,f)
!else if(Iflag_scheme.eq.100) then   !中心差分格式必须同weno混合才能使用
!    call center6(Ka,Kb,fpa,fma,f)
else if(Iflag_scheme.eq.11) then
    call Hweno5_c6(Ka,Kb,fpa,fma,f)
endif


return
end