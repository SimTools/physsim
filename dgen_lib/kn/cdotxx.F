C
C **********************************************************************
C
      SUBROUTINE CDOTxx(V1,V2 , SCALAR)
C
C This subroutine computes the scalar product of two COMPLEX*16 vectors.
C For the real vectors (momenta), use RDOTxx subroutine.
C
C INPUT:
C       complex V1(0:3)        : vector
C       complex V2(0:3)        : vector
C
C OUTPUT:
C       complex SCALAR         : scalar product of two vectors
C
      COMPLEX*16 V1(0:3),V2(0:3),SCALAR
C
      SCALAR =  V1(0)*V2(0)-V1(1)*V2(1)
     &         -V1(2)*V2(2)-V1(3)*V2(3)
      RETURN
      END
