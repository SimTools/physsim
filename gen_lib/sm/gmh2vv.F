C* (Inputs)
C*    AMH    : (R*4) : Higss mass.
C*    AMV    : (R*4) : vector boson mass.
C*    AMZ    : (R*4) : Z mass.
C* (Output)
C*    GM     : (R*4) : (1/1)*( Higss partical width ) for V = Z.
C*                   : (1/2)*( Higss partical width ) for V = W.
C*
 
      SUBROUTINE GMH2VV(AMH,AMV,AMZ,GM)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
#include "../sm/include/SMCONS.inc"
      REAL*4 AMH, AMV, AMZ, GM
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      GF     = xPI*xALF/xSIN2W/xCOS2W/xSQRT2/AMZ**2
      FACT   = GF/16/xSQRT2/xPI
C--
C  Calculate width.
C--
      EW = AMH/2
      BT = (EW-AMV)*(EW+AMV)
      IF ( BT.LE.0. ) THEN
         GM = 0
         RETURN
      ELSE
         BT   = SQRT(BT)/EW
      ENDIF
      GM = FACT*AMH**3*BT*( 1 - 4*(AMV/AMH)**2 + 12*(AMV/AMH)**4 )
C--
C  That's it.
C--
      RETURN
      END
