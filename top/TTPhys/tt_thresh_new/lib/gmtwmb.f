CC**********************************************************************
C*  real*8  function GMTWMB( GF, mt, mW, mb )
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION GMTWMB( GF, MT, MW, MB )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       GF, MT, MW, MB
C--
      REAL   *8       K, EB, PI, RB, RW
C
C========< Entry Point >================================================
C
C--
C  Calculate the lowest Gamma_t with mb.
C--
      PI = ACOS(-1.D0)
      RB = MB**2/MT**2
      RW = MW**2/MT**2
      EB = ( MT**2 - MW**2 + MB**2 )/2/MT
      K  = SQRT( EB**2 - MB**2 )
      GMTWMB = GF*MT**3/8/SQRT(2.D0)/PI * 2*K/MT
     .            *( (1-RB)**2 + (1+RB)*RW - 2*RW**2 )
C--
C  That's it.
C--
      RETURN
      END
