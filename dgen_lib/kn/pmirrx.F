C
C ----------------------------------------------------------------------
C
      SUBROUTINE PMIRRx(P , PMIRR)
C
C This subroutine computes the real vector whose signs of the spacial
C components are reversed from the ones of the original vector.
C
C INPUT:
C       real    P(0:3)         : real vector
C
C OUTPUT:
C       real    PMIRR(0:3)     : mirror vector
C
      REAL*8    P(0:3),PMIRR(0:3)
C
      PMIRR(0)= P(0)
      PMIRR(1)=-P(1)
      PMIRR(2)=-P(2)
      PMIRR(3)=-P(3)
      RETURN
      END
