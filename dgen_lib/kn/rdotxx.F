C
C ----------------------------------------------------------------------
C
      SUBROUTINE RDOTxx(P1,P2 , SCALAR)
C
C This subroutine computes the scalar product of two REAL*8 vectors.
C To get an inner-product of two momenta, call this subroutine.
C For the complex vectors, use CDOTxx subroutine.
C
C INPUT:
C       real    P1(0:3)        : vector
C       real    P2(0:3)        : vector
C
C OUTPUT:
C       real    SCALAR         : scalar product of two vectors
C
      REAL*8 P1(0:3),P2(0:3),SCALAR
C
      SCALAR =  P1(0)*P2(0)-P1(1)*P2(1)
     &         -P1(2)*P2(2)-P1(3)*P2(3)
      RETURN
      END
