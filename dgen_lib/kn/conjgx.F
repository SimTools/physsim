C
C ======================================================================
C
      SUBROUTINE CONJGx(V , VCONJG)
C
C This subroutine computes the conjugation of a complex vector.
C
C INPUT:
C       complex V(0:3)         : complex vector
C
C OUTPUT:
C       complex VCONJG(0:3)    : conjugation V*
C
      COMPLEX*16 V(0:3),VCONJG(0:3)
C
      VCONJG(0)=CONJG(V(0))
      VCONJG(1)=CONJG(V(1))
      VCONJG(2)=CONJG(V(2))
      VCONJG(3)=CONJG(V(3))
      RETURN
      END
