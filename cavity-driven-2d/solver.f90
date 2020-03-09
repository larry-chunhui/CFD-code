subroutine solver(I,J,b1,b2,b3,b4,b5)
include 'incl.h'
real*8 b1,b2,b3,b4,b5
integer I,J,k
!
select case(Id_Escheme)
case(1)  !1 upwind∏Ò Ω
	call one_upwind
case(3)
	call quick_upwind
end select

selectcase(Id_Dscheme)
case(2) 
	call two_diffuse
case(4)
	call four_diffuse
end select
!***************************************************
do k=0,1
      AEU(k)=AEU_con(k)+AE_diff
      AWU(k)=AWU_con(k)+AW_diff
      ANU(k)=ANU_con(k)+AN_diff
      ASU(k)=ASU_con(k)+AS_diff
 
      AEV(k)=AEV_con(k)+AE_diff
      AWV(k)=AWV_con(k)+AW_diff
      ANV(k)=ANV_con(k)+AN_diff
      ASV(k)=ASV_con(k)+AS_diff

      AEEU(k)=AEEU_con(k)+AEE_diff
      AWWU(k)=AWWU_con(k)+AWW_diff
      ANNU(k)=ANNU_con(k)+ANN_diff
      ASSU(k)=ASSU_con(k)+ASS_diff
 
      AEEV(k)=AEEV_con(k)+AEE_diff
      AWWV(k)=AWWV_con(k)+AWW_diff
      ANNV(k)=ANNV_con(k)+ANN_diff
      ASSV(k)=ASSV_con(k)+ASS_diff

      APU(k)=AEU(k)+AWU(k)+ANU(k)+ASU(k)+AEEU(k)+AWWU(k)+ANNU(k)+ASSU(k)
      APV(k)=AEV(k)+AWV(k)+ANV(k)+ASV(k)+AEEV(k)+AWWV(k)+ANNV(k)+ASSV(k)
enddo
!***************************************************
       BU(1)=( P(I,J)-P(I+1,J) )*DY  
       BU(0)=( P(I-1,J)-P(I,J) )*DY  
       BV(1)=( P(I,J)-P(I,J+1) )*DX 
       BV(0)=( P(I,J-1)-P(I,J) )*DX 
!U(i-1,j)
      b1=-APU(0)*U(I-1,J)+AEU( 0)*U(I,  J)+AWU( 0)*U(I-2,J)+ANU( 0)*U(I-1,J+1)+ASU( 0)*U(I-1,J-1)+BU(0) &
						 +AEEU(0)*U(I+1,J)+AWWU(0)*U(I-3,J)+ANNU(0)*U(I-1,J+2)+ASSU(0)*U(I-1,J-2)
!U(i,j)   
      b2=-APU(1)*U(I,  J)+AEU( 1)*U(I+1,J) +AWU(1)*U(I-1,J) +ANU(1)*U(I,J+1)+ASU( 1)*U(I,J-1)+BU(1) &
						 +AEEU(1)*U(I+2,J)+AWWU(1)*U(I-2,J)+ANNU(1)*U(I,J+2)+ASSU(1)*U(I,J-2)
!V(i,j-1)
      b3=-APV(0)*V(I,J-1)+AEV( 0)*V(I+1,J-1)+AWV( 0)*V(I-1,J-1)+ANV( 0)*V(I,J  )+ASV( 0)*V(I,J-2)+BV(0) &
						 +AEEV(0)*V(I+2,J-1)+AWWV(0)*V(I-2,J-1)+ANNV(0)*V(I,J+1)+ASSV(0)*V(I,J-3)
!V(i,j)
      b4=-APV(1)*V(I,  J)+AEV( 1)*V(I+1,J)+AWV( 1)*V(I-1,J)+ANV( 1)*V(I,  J+1)+ASV( 1)*V(I,  J-1)+BV(1) &
						 +AEEV(1)*V(I+2,J)+AWWV(1)*V(I-2,J)+ANNV(1)*V(I,  J+2)+ASSV(1)*V(I,  J-2)
      b5=( U(I-1,J)-U(I,J) )*DY + ( V(I,J-1)-V(I,J) )*DX
!a still problem!!!!***************************8

return
end