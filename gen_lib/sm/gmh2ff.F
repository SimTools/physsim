C* (Inputs)
C*    AMH    : (R*4) : Higss mass.
C*    AMF    : (R*4) : fermion mass.
C*    AMZ    : (R*4) : Z mass.
C* (Output)
C*    GM     : (R*4) : (1/1)*( Higss partical width ) for F = L.
C*                   : (1/3)*( Higss partical width ) for F = Q.
C*
 
      SUBROUTINE GMH2FF(AMH,AMF,AMZ,GM)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
#include "../sm/include/SMCONS.inc"
      REAL*4 AMH, AMF, AMZ, GM
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      GF     = xPI*xALF/xSIN2W/xCOS2W/xSQRT2/AMZ**2
      FACT   = GF/4/xSQRT2/xPI
C--
C  Calculate width.
C--
      EF = AMH/2
      BT = (EF-AMF)*(EF+AMF)
      IF ( BT.LE.0. ) THEN
         GM = 0
         RETURN
      ELSE
         BT   = SQRT(BT)/EF
      ENDIF
      GM = FACT*AMF*AMF*AMH*BT**3
C--
C  That's it.
C--
      RETURN
      END
