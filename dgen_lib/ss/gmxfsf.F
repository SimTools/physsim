      SUBROUTINE GMXFSF(GCFSF,AMX,AMF,AMSF,GM)

      IMPLICIT REAL*8 ( A-H, O-Z )
#include "../sm/include/SMCONS.inc"
      COMPLEX*16 GCFSF(2)
      REAL* 8    AMX, AMF, AMSF, GM
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT( 1 - 2*(X1+X2) + (X1-X2)**2 )
C
C========< Entry Point >================================================
C
C--
C  Calculate width.
C--
      IF ( AMF+AMSF.GE.AMX ) THEN
         GM = 0
         RETURN
      ENDIF
      X1  = (AMF/AMX)**2
      X2  = (AMSF/AMX)**2
      FAC = BETA(X1,X2)/(32*xPI*AMX)
      GM  = FAC*( (ABS(GCFSF(1))**2+ABS(GCFSF(2))**2)
     .		 *( (AMX-AMSF)*(AMX+AMSF) + AMF*AMF )
     .		+ 4*ABS(GCFSF(1)*GCFSF(2))*AMX*AMF )
C--
C  That's it.
C--
      RETURN
      END
