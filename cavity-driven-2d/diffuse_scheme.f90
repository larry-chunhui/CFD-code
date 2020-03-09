
subroutine two_diffuse
include 'incl.h'

AP_diff=4.0D0/RE
AE_diff=1.0d0/RE
AW_diff=1.0d0/RE
AN_diff=1.0d0/RE
AS_diff=1.0d0/RE

AEE_diff=0.0D0
AWW_diff=0.0D0
ANN_diff=0.0D0
ASS_diff=0.0D0
return
end



subroutine four_diffuse
include 'incl.h'
AP_diff=-5.0D0/RE

AW_diff=4.0d0/3.0d0/RE
AE_diff=4.0d0/3.0d0/RE
AS_diff=4.0d0/3.0d0/RE
AN_diff=4.0d0/3.0d0/RE

AWW_diff=-1.0D0/12.0D0/RE
AEE_diff=-1.0D0/12.0D0/RE
ASS_diff=-1.0D0/12.0D0/RE
ANN_diff=-1.0D0/12.0D0/RE
return
end