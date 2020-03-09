subroutine H_smooth(a0,a1,a2,tao)
    use var
    implicit none
    real*8 a0,a1,a2,tao
    real*8 ep
    
    ep=1.d-6
    if(Iflag_smooth.eq.0) then
        tao=(dabs(a0-a2)/(MIN(a0,a1,a2)+ep))**2
    else if(Iflag_smooth.eq.1) then
        tao=((max(a0,a1,a2)-min(a0,a1,a2))/(min(a0,a1,a2)+ep))**2
    endif
    
    
    return
    end
    