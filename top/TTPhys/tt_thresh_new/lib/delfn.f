CC*********************************************************************
C*  REAL*8 FUNCTION  DELFN(X)
CC*********************************************************************
 
      DOUBLE PRECISION FUNCTION DELFN(X)
 
      REAL*8   X
C
C===========< Entry Point >=============================================
C
C--
C  APPROXIMATED DELTA FUNCTION.
C--
      DELFN = (SIN(X)-X*COS(X))/X
C--
C  That's it.
C--
      RETURN
      END
