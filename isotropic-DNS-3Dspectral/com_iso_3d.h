
	IMPLICIT NONE

	INTEGER  E_ijk(3,2)
	COMMON/E_ijk/E_ijk
        
    REAL*8 XL,YL,ZL,Re,DT
	COMMON /GP1/XL,YL,ZL,Re,DT

	INTEGER ICONT,NSTEPS,IDUMFRQ,IPRNFRQ,IWRITFRQ
	COMMON /GP2/ICONT,NSTEPS,IDUMFRQ,IPRNFRQ,IWRITFRQ
                          

	INTEGER IFORCE,KFMAX
	COMMON /GP3/ IFORCE,KFMAX


