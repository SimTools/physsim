C* (Inputs)
C*    AMHM   : (R*4) : H- mass.
C*    AMU    : (R*4) : U-type fermion mass.
C*    AMD    : (R*4) : D-type fermion mass.
C*    AMZ    : (R*4) : Z mass.
C*    TNB    : (R*4) : tan(beta).
C* (Output)
C*    GM     : (R*4) : (1/1)*( Higss partical width ) for leptons.
C*                   : (1/3)*( Higss partical width ) for quarks.
C*
 
      SUBROUTINE GMH2UD(AMHM,AMU,AMD,AMZ,TNB,GM)
 
      IMPLICIT REAL*8 ( A-H, O-Z    )
#include "../sm/include/SMCONS.inc"
      REAL*8 AMHM, AMU, AMD, AMZ, TNB, GM
C
      BETA(X1,X2) = SQRT( 1 - 2*(X1+X2) + (X1-X2)**2 )
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      FACT   = xALF/xSIN2W/xCOS2W/AMZ**2/16
C--
C  Calculate width.
C--
      IF ( AMHM.LE.AMU+AMD) THEN
         GM = 0
         RETURN
      ELSE
         BT = BETA((AMU/AMHM)**2,(AMD/AMHM)**2)
      ENDIF
      A2 = ( AMU/TNB + AMD*TNB )**2
      B2 = ( AMU/TNB - AMD*TNB )**2
      GM = FACT*( A2*(AMHM-AMU-AMD)*(AMHM+AMU+AMD)
     .          + B2*(AMHM-AMU+AMD)*(AMHM+AMU-AMD) )*BT/AMHM
C--
C  That's it.
C--
      RETURN
      END
