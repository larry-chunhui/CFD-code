subroutine one_upwind(I,J,b1,b2,b3,b4,b5)
include 'incl.h'
real*8 b1,b2,b3,b4,b5
integer I,J
      AEU(1)=1.0d0/RE+max(  0.0d0,-fe_u(1) )
      AWU(1)=1.0d0/RE+max(  0.0d0, fw_u(1) )
      ANU(1)=1.0d0/RE+max(  0.0d0,-fn_u(1) )
      ASU(1)=1.0d0/RE+max(  0.0d0, fs_u(1) )
      APU(1)=AEU(1)+AWU(1)+ANU(1)+ASU(1)
       BU(1)=( P(I,J)-P(I+1,J) )*DY    
!c for u(i-1,j)
      AEU(0)=1.0d0/RE+max(  0.0d0,-fe_u(0) )
      AWU(0)=1.0d0/RE+max(  0.0d0, fw_u(0) )
      ANU(0)=1.0d0/RE+max(  0.0d0,-fn_u(0) )
      ASU(0)=1.0d0/RE+max(  0.0d0, fs_u(0) )
      APU(0)=AEU(0)+AWU(0)+ANU(0)+ASU(0)
       BU(0)=( P(I-1,J)-P(I,J) )*DY  
!c for v(i,j  )
      AEV(1)=1.0d0/RE+max(  0.0d0,-fe_v(1) )
      AWV(1)=1.0d0/RE+max(  0.0d0, fw_v(1) )
      ANV(1)=1.0d0/RE+max(  0.0d0,-fn_v(1) )
      ASV(1)=1.0d0/RE+max(  0.0d0, fs_v(1) )
      APV(1)=AEV(1)+AWV(1)+ANV(1)+ASV(1)
       BV(1)=( P(I,J)-P(I,J+1) )*DX 	
!c for v(i,j-1)
      AEV(0)=1.0d0/RE+max(  0.0d0,-fe_v(0)  )
      AWV(0)=1.0d0/RE+max(  0.0d0, fw_v(0)  )
      ANV(0)=1.0d0/RE+max(  0.0d0,-fn_v(0)  )
      ASV(0)=1.0d0/RE+max(  0.0d0, fs_v(0)  )
      APV(0)=AEV(0)+AWV(0)+ANV(0)+ASV(0)
       BV(0)=( P(I,J-1)-P(I,J) )*DX 
!U(i-1,j)
      b1=-APU(0)*U(I-1,J)+AEU(0)*U(I,  J)+AWU(0)*U(I-2,J)+ANU(0)*U(I-1,J+1)+ASU(0)*U(I-1,J-1)+BU(0) 
!U(i,j)   
      b2=-APU(1)*U(I,  J)+AEU(1)*U(I+1,J)+AWU(1)*U(I-1,J)+ANU(1)*U(I,  J+1)+ASU(1)*U(I,  J-1)+BU(1)
!V(i,j-1)
      b3=-APV(0)*V(I,J-1)+AEV(0)*V(I+1,J-1)+AWV(0)*V(I-1,J-1)+ANV(0)*V(I,J)+ASV(0)*V(I,J-2)+BV(0)
!V(i,j)
      b4=-APV(1)*V(I,  J)+AEV(1)*V(I+1,J)+AWV(1)*V(I-1,J)+ANV(1)*V(I,  J+1)+ASV(1)*V(I,  J-1)+BV(1)
      b5=( U(I-1,J)-U(I,J) )*DY + ( V(I,J-1)-V(I,J) )*DX

return
end