C
C ----------------------------------------------------------------------
C
      SUBROUTINE PSUMxx(P1,P2 , PSUM)
C
C This subroutine computes the sum of two real vectors.
C
C INPUT:
C       real    P1(0:3)        : real vector
C       real    P2(0:3)        : real vector
C
C OUTPUT:
C       real    PSUM(0:3)      : summation P1+P2
C
      REAL*8    P1(0:3),P2(0:3),PSUM(0:3)
C
      PSUM(0)=P1(0)+P2(0)
      PSUM(1)=P1(1)+P2(1)
      PSUM(2)=P1(2)+P2(2)
      PSUM(3)=P1(3)+P2(3)
      RETURN
      END
