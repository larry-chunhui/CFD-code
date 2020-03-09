subroutine hybrid(I,J,b1,b2,b3,b4,b5)
include 'incl.h'
real*8 b1,b2,b3,b4,b5
integer I,J
integer k
!write(*,*)"------------------------------"
!write(*,*)max(100.0d0,2.0d0,8.0d0)
!stop
do k=0,1
      AEU(k)=max( -fe_u(k), 1.0d0/RE-fe_u(k)/2.0d0, 0.0d0 )
      AWU(k)=max(  fw_u(k), 1.0d0/RE+fw_u(k)/2.0d0, 0.0d0 )
      ANU(k)=max( -fn_u(k), 1.0d0/RE-fn_u(k)/2.0d0, 0.0d0 )
      ASU(k)=max(  fs_u(k), 1.0d0/RE+fs_u(k)/2.0d0, 0.0d0 )
      APU(k)=AEU(k)+AWU(k)+ANU(k)+ASU(k)

      AEV(k)=max( -fe_v(k), 1.0d0/RE-fe_v(k)/2.0d0, 0.0d0 )
      AWV(k)=max(  fw_v(k), 1.0d0/RE+fw_v(k)/2.0d0, 0.0d0 )
      ANV(k)=max( -fn_v(k), 1.0d0/RE-fn_v(k)/2.0d0, 0.0d0 )
      ASV(k)=max(  fs_v(k), 1.0d0/RE+fs_v(k)/2.0d0, 0.0d0 )
      APV(k)=AEV(k)+AWV(k)+ANV(k)+ASV(k)
enddo
       BU(1)=( P(I,J)-P(I+1,J) )*DY    
!c for u(i-1,j)

       BU(0)=( P(I-1,J)-P(I,J) )*DY  
!c for v(i,j  )

       BV(1)=( P(I,J)-P(I,J+1) )*DX 	
!c for v(i,j-1)

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