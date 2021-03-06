C     THIS FUNCTION SUPPLIES THE SQUARE OF THE MATRIX ELEMENTS
C     AND IS USED IN COMBINATION WITH THE BASES MONTECARLO ROUTINE
C     THE FUNCTION IS CALLED FROM FUNC OF THE BASES ROUTINE
C                                             20.DEZ.,  1986 M.KURODA
               SUBROUTINE FUNC3(TAU4,TAU5)
 
 
 
           IMPLICIT REAL*8(A-H,O-Z)
 
#include "KINEM.inc"
#include "KINEM2.inc"
 
 
 
      ANS4=-P1Q2*P3Q1**2-P1Q2*P3Q1*EME2-3.*P1Q2*P3Q1*UME2+2.*
     . P3Q2**2*EME2
      ANS3=-4.*P1P3**3+4.*P1P3**2*P1Q2+4.*P1P3**2*P3Q2-4.*P1P3**
     . 2*P3Q1-2.*P1P3*S2P1*S2Q2-2.*P1P3*S2P1*S2Q1+4.*P1P3*S2P3**
     . 2+2.*P1P3*S2P3*S2Q2+8.*P1P3*S2P3*S2Q1-2.*P1P3*S2Q2**2-4.*
     . P1P3*S2Q2*S2Q1-2.*P1P3*S2Q1**2-2.*P1P3*P1Q2**2-8.*P1P3*
     . P1Q2*P3Q2-4.*P1P3*P1Q2*UME2+4.*P1P3*P3Q2*P3Q1-4.*P1P3*
     . P3Q2*EME2-2.*P1P3*P3Q1*EME2+2.*S2P1**2*P3Q2+2.*S2P1**2*
     . P3Q1+4.*S2P1*S2P3*P1Q2+8.*S2P1*S2P3*P3Q2+4.*S2P1*S2P3*
     . EME2+4.*S2P1*S2P3*UME2+2.*S2P1*S2Q2*P3Q2-2.*S2P1*S2Q2*
     . UME2+4.*S2P1*S2Q1*P3Q2+2.*S2P1*S2Q1*P3Q1-16.*S2P3**2*P1Q2
     . +8.*S2P3*S2Q2*P1Q2-2.*S2P3*S2Q2*P3Q2+2.*S2P3*S2Q2*P3Q1+2.
     . *S2P3*S2Q2*UME2+2.*S2P3*S2Q1*P1Q2+2.*S2P3*S2Q1*P3Q2+4.*
     . S2P3*S2Q1*P3Q1+2.*S2P3*S2Q1*EME2+4.*S2P3*S2Q1*UME2-2.*
     . S2Q2**2*P3Q1-S2Q2**2*UME2+2.*S2Q2*S2Q1*P3Q2-2.*S2Q2*S2Q1*
     . P3Q1-2.*S2Q2*S2Q1*UME2+2.*S2Q1**2*P3Q2-8.*S22*P1Q2*P3Q2-
     . 5.*S22*P1Q2*P3Q1+S22*P3Q2**2+6.*P1Q2**2*P3Q2+2.*P1Q2**2*
     . P3Q1+6.*P1Q2**2*UME2-4.*P1Q2*P3Q2*P3Q1-2.*P1Q2*P3Q2*EME2-
     . 2.*P1Q2*P3Q2*UME2+ANS4
      ANS2=4.*Q12*Q22*ANS3
      ANS8=-4.*S2Q2**2*P1Q2*P3Q1-2.*S2Q2**2*P1Q2*UME2+S2Q2**2*
     . P3Q1*EME2+4.*S2Q2*S2Q1*P1Q2*P3Q2-S2Q2*S2Q1*P1Q2*UME2-S2Q2
     . *S2Q1*P3Q2*EME2-6.*S22*P1Q2**2*P3Q2-2.*S22*P1Q2**2*P3Q1+
     . 2.*S22*P1Q2*P3Q2**2+S22*P1Q2*P3Q2*P3Q1
      ANS7=8.*P1P3**2*S2P1*S2P3-4.*P1P3**2*S2P1*S2Q2-4.*P1P3**2*
     . S2P3*S2Q2+2.*P1P3**2*S2Q2**2+4.*P1P3*S2P1**2*P3Q2-8.*P1P3
     . *S2P1*S2P3*P1Q2-4.*P1P3*S2P1*S2P3*P3Q2+8.*P1P3*S2P1*S2P3*
     . P3Q1-4.*P1P3*S2P1*S2Q2*P3Q1+2.*P1P3*S2P1*S2Q1*P3Q2+8.*
     . P1P3*S2P3**2*P1Q2+4.*P1P3*S2P3*S2Q2*P1Q2-4.*P1P3*S2P3*
     . S2Q2*P3Q1+4.*P1P3*S2P3*S2Q1*P1Q2-4.*P1P3*S2Q2**2*P1Q2+2.*
     . P1P3*S2Q2**2*P3Q1-4.*P1P3*S2Q2*S2Q1*P1Q2+2.*P1P3*S22*P1Q2
     . **2-2.*S2P1**2*P3Q2**2+2.*S2P1**2*P3Q2*P3Q1+12.*S2P1*S2P3
     . *P1Q2*P3Q2-4.*S2P1*S2P3*P1Q2*P3Q1-4.*S2P1*S2P3*P3Q2*P3Q1+
     . 4.*S2P1*S2P3*P3Q2*EME2+4.*S2P1*S2P3*P3Q1*EME2+4.*S2P1*S2Q2
     . *P1Q2*P3Q2+4.*S2P1*S2Q2*P1Q2*P3Q1-2.*S2P1*S2Q2*P3Q1*EME2-
     . 2.*S2P1*S2Q1*P3Q2**2+2.*S2P1*S2Q1*P3Q2*EME2-12.*S2P3**2*
     . P1Q2**2+4.*S2P3**2*P1Q2*P3Q1+6.*S2P3*S2Q2*P1Q2**2-4.*S2P3
     . *S2Q2*P1Q2*P3Q2+2.*S2P3*S2Q2*P1Q2*P3Q1+2.*S2P3*S2Q2*P1Q2*
     . EME2+4.*S2P3*S2Q2*P1Q2*UME2-2.*S2P3*S2Q2*P3Q2*EME2-2.*
     . S2P3*S2Q2*P3Q1*EME2+4.*S2P3*S2Q1*P1Q2*P3Q2+2.*S2P3*S2Q1*
     . P1Q2*P3Q1+2.*S2P3*S2Q1*P1Q2*EME2+2.*S2P3*S2Q1*P1Q2*UME2+
     . ANS8
      ANS6=8.*Q12*ANS7
      ANS14=-4.*S2P1*S2P3*P3Q1*EME2+4.*S2P1*S2Q2*P1Q2*P3Q2+2.*
     . S2P1*S2Q2*P1Q2*P3Q1+2.*S2P1*S2Q2*P3Q1*EME2+4.*S2P1*S2Q1*
     . P1Q2*P3Q2+2.*S2P1*S2Q1*P1Q2*P3Q1-2.*S2P1*S2Q1*P3Q2**2-2.*
     . S2P1*S2Q1*P3Q2*P3Q1-2.*S2P1*S2Q1*P3Q2*EME2-12.*S2P3**2*
     . P1Q2**2+6.*S2P3*S2Q2*P1Q2**2+4.*S2P3*S2Q2*P1Q2*P3Q1-2.*
     . S2P3*S2Q2*P1Q2*EME2-2.*S2P3*S2Q2*P3Q2*EME2-2.*S2P3*S2Q2*
     . P3Q1*EME2+2.*S2P3*S2Q1*P1Q2**2+4.*S2P3*S2Q1*P1Q2*P3Q2+6.*
     . S2P3*S2Q1*P1Q2*P3Q1-2.*S2P3*S2Q1*P1Q2*EME2+2.*S2P3*S2Q1*
     . P1Q2*UME2-4.*S2P3*S2Q1*P3Q2*EME2-4.*S2P3*S2Q1*P3Q1*EME2-
     . 2.*S2Q2**2*P1Q2*P3Q1+S2Q2**2*P3Q1*EME2+2.*S2Q2*S2Q1*P1Q2*
     . P3Q2-2.*S2Q2*S2Q1*P1Q2*P3Q1-S2Q2*S2Q1*P1Q2*UME2-S2Q2*S2Q1
     . *P3Q2*EME2+2.*S2Q2*S2Q1*P3Q1*EME2+2.*S2Q1**2*P1Q2*P3Q2-2.
     . *S2Q1**2*P3Q2*EME2-6.*S22*P1Q2**2*P3Q2-4.*S22*P1Q2**2*
     . P3Q1-S22*P1Q2*P3Q2*P3Q1+4.*P1Q2**3*P3Q2+2.*P1Q2**3*P3Q1+
     . 4.*P1Q2**3*UME2-4.*P1Q2**2*P3Q2*P3Q1-2.*P1Q2**2*P3Q1**2-2.
     . *P1Q2**2*P3Q1*UME2+4.*P1Q2*P3Q2**2*EME2+4.*P1Q2*P3Q2*P3Q1
     . *EME2+2.*P1Q2*P3Q1**2*EME2
      ANS13=-8.*P1P3**2*S2P1*S2P3+4.*P1P3**2*S2P1*S2Q2-4.*P1P3**
     . 2*S2P3*S2Q2-8.*P1P3**2*S2P3*S2Q1+2.*P1P3**2*S2Q2**2+4.*
     . P1P3**2*S2Q2*S2Q1+8.*P1P3**2*P1Q2*P3Q2+4.*P1P3**2*P1Q2*
     . P3Q1-4.*P1P3*S2P1**2*P3Q2+8.*P1P3*S2P1*S2P3*P1Q2-4.*P1P3*
     . S2P1*S2P3*P3Q2-8.*P1P3*S2P1*S2P3*P3Q1-4.*P1P3*S2P1*S2Q2*
     . P1Q2+4.*P1P3*S2P1*S2Q2*P3Q1-4.*P1P3*S2P1*S2Q1*P1Q2-6.*
     . P1P3*S2P1*S2Q1*P3Q2+8.*P1P3*S2P3**2*P1Q2+4.*P1P3*S2P3*
     . S2Q2*P1Q2-4.*P1P3*S2P3*S2Q2*P3Q1+12.*P1P3*S2P3*S2Q1*P1Q2-
     . 4.*P1P3*S2P3*S2Q1*P3Q2-8.*P1P3*S2P3*S2Q1*P3Q1-4.*P1P3*S2Q2
     . **2*P1Q2+2.*P1P3*S2Q2**2*P3Q1-6.*P1P3*S2Q2*S2Q1*P1Q2+4.*
     . P1P3*S2Q2*S2Q1*P3Q1-2.*P1P3*S2Q1**2*P1Q2-2.*P1P3*S2Q1**2*
     . P3Q2-2.*P1P3*S22*P1Q2**2-8.*P1P3*P1Q2**2*P3Q2-4.*P1P3*
     . P1Q2**2*P3Q1-4.*P1P3*P1Q2**2*UME2+8.*P1P3*P1Q2*P3Q2*P3Q1+
     . 4.*P1P3*P1Q2*P3Q1**2+4.*S2P1**2*P1Q2*P3Q2+4.*S2P1**2*P1Q2*
     . P3Q1-2.*S2P1**2*P3Q2**2-2.*S2P1**2*P3Q2*P3Q1+4.*S2P1*S2P3
     . *P1Q2**2+12.*S2P1*S2P3*P1Q2*P3Q2+4.*S2P1*S2P3*P1Q2*P3Q1-
     . 4.*S2P1*S2P3*P3Q2*EME2+ANS14
      ANS12=8.*Q22*ANS13
      ANS18=-S2Q2*S2Q1*P3Q2*EME2-2.*S2Q1**2*P1Q2*P3Q2-6.*S22*
     . P1Q2**2*P3Q2-2.*S22*P1Q2**2*P3Q1+S22*P1Q2*P3Q2*P3Q1
      ANS17=8.*P1P3**2*S2P1*S2P3-4.*P1P3**2*S2P1*S2Q2-4.*P1P3**2
     . *S2P3*S2Q2+2.*P1P3**2*S2Q2**2+4.*P1P3*S2P1**2*P3Q2-8.*
     . P1P3*S2P1*S2P3*P1Q2-4.*P1P3*S2P1*S2P3*P3Q2+8.*P1P3*S2P1*
     . S2P3*P3Q1-4.*P1P3*S2P1*S2Q2*P3Q1+2.*P1P3*S2P1*S2Q1*P3Q2+
     . 8.*P1P3*S2P3**2*P1Q2+4.*P1P3*S2P3*S2Q2*P1Q2-4.*P1P3*S2P3*
     . S2Q2*P3Q1+4.*P1P3*S2P3*S2Q1*P1Q2-4.*P1P3*S2Q2**2*P1Q2+2.*
     . P1P3*S2Q2**2*P3Q1-4.*P1P3*S2Q2*S2Q1*P1Q2+2.*P1P3*S22*P1Q2
     . **2-2.*S2P1**2*P3Q2**2+2.*S2P1**2*P3Q2*P3Q1+12.*S2P1*S2P3
     . *P1Q2*P3Q2-4.*S2P1*S2P3*P1Q2*P3Q1-4.*S2P1*S2P3*P3Q2*P3Q1+
     . 4.*S2P1*S2P3*P3Q2*EME2+4.*S2P1*S2P3*P3Q1*EME2+4.*S2P1*S2Q2
     . *P1Q2*P3Q2+4.*S2P1*S2Q2*P1Q2*P3Q1-2.*S2P1*S2Q2*P3Q1*EME2-
     . 2.*S2P1*S2Q1*P3Q2**2+2.*S2P1*S2Q1*P3Q2*EME2-12.*S2P3**2*
     . P1Q2**2+4.*S2P3**2*P1Q2*P3Q1+6.*S2P3*S2Q2*P1Q2**2-2.*S2P3
     . *S2Q2*P1Q2*P3Q1+2.*S2P3*S2Q2*P1Q2*EME2-2.*S2P3*S2Q2*P3Q2*
     . EME2-2.*S2P3*S2Q2*P3Q1*EME2+2.*S2P3*S2Q1*P1Q2*P3Q2-2.*
     . S2P3*S2Q1*P1Q2*P3Q1+2.*S2P3*S2Q1*P1Q2*EME2-2.*S2P3*S2Q1*
     . P1Q2*UME2+S2Q2**2*P3Q1*EME2+2.*S2Q2*S2Q1*P1Q2*P3Q1+S2Q2*
     . S2Q1*P1Q2*UME2+ANS18
      ANS16=8.*Q32*ANS17
      ANS19=32.*P1Q2*(-4.*P1P3**2*S2P3*S2Q2-4.*P1P3**2*S2P3*S2Q1
     . +2.*P1P3**2*S2Q2**2+2.*P1P3**2*S2Q2*S2Q1-4.*P1P3*S2P1*
     . S2P3*P3Q2-2.*P1P3*S2P1*S2Q1*P3Q2+4.*P1P3*S2P3**2*P1Q2+2.*
     . P1P3*S2P3*S2Q2*P1Q2-4.*P1P3*S2P3*S2Q2*P3Q1+4.*P1P3*S2P3*
     . S2Q1*P1Q2-2.*P1P3*S2P3*S2Q1*P3Q2-4.*P1P3*S2P3*S2Q1*P3Q1-
     . 2.*P1P3*S2Q2**2*P1Q2+2.*P1P3*S2Q2**2*P3Q1-2.*P1P3*S2Q2*
     . S2Q1*P1Q2+2.*P1P3*S2Q2*S2Q1*P3Q1-P1P3*S2Q1**2*P3Q2-2.*
     . S2P1**2*P3Q2**2+6.*S2P1*S2P3*P1Q2*P3Q2-2.*S2P1*S2P3*P3Q2*
     . P3Q1+2.*S2P1*S2Q2*P1Q2*P3Q2+2.*S2P1*S2Q2*P1Q2*P3Q1-2.*
     . S2P1*S2Q1*P3Q2**2-S2P1*S2Q1*P3Q2*P3Q1-4.*S2P3**2*P1Q2**2+
     . 2.*S2P3**2*P1Q2*P3Q1+2.*S2P3*S2Q2*P1Q2**2+S2P3*S2Q2*P1Q2*
     . P3Q1-2.*S2P3*S2Q2*P3Q2*EME2-2.*S2P3*S2Q2*P3Q1*EME2+3.*
     . S2P3*S2Q1*P1Q2*P3Q2+2.*S2P3*S2Q1*P1Q2*P3Q1-2.*S2P3*S2Q1*
     . P3Q2*EME2-2.*S2P3*S2Q1*P3Q1*EME2-S2Q2**2*P1Q2*P3Q1+S2Q2**
     . 2*P3Q1*EME2+S2Q2*S2Q1*P1Q2*P3Q2-S2Q2*S2Q1*P3Q2*EME2+S2Q2*
     . S2Q1*P3Q1*EME2-S2Q1**2*P3Q2*EME2-2.*S22*P1Q2**2*P3Q2-S22*
     . P1Q2**2*P3Q1)
      ANS15=2.*Q32**3*(P1P3*S22-2.*S2P1*S2P3)+4.*Q32**2*(-4.*
     . P1P3*S2P1*S2P3+4.*P1P3*S2P3*S2Q2-2.*P1P3*S2Q2**2-2.*P1P3*
     . S2Q2*S2Q1+2.*P1P3*S22*P1Q2-P1P3*S22*P3Q2+4.*S2P1*S2P3*
     . P3Q2+2.*S2P1*S2P3*P3Q1-2.*S2P1*S2P3*EME2+2.*S2P1*S2P3*
     . UME2+2.*S2P1*S2Q2*P3Q2-S2P1*S2Q2*UME2+2.*S2P1*S2Q1*P3Q2-
     . 4.*S2P3**2*P1Q2+2.*S2P3*S2Q2*P1Q2+S2P3*S2Q2*EME2+2.*S2P3*
     . S2Q1*P1Q2-2.*S22*P1Q2*P3Q2-S22*P1Q2*P3Q1)+ANS16+ANS19
      ANS11=4.*Q22*Q32*(-4.*P1P3**3+4.*P1P3**2*P1Q2+4.*P1P3**2*
     . P3Q2-4.*P1P3**2*P3Q1-2.*P1P3*S2P1*S2Q2-2.*P1P3*S2P1*S2Q1+
     . 4.*P1P3*S2P3**2+2.*P1P3*S2P3*S2Q2+8.*P1P3*S2P3*S2Q1-2.*
     . P1P3*S2Q2**2-4.*P1P3*S2Q2*S2Q1-2.*P1P3*S2Q1**2-2.*P1P3*
     . P1Q2**2-8.*P1P3*P1Q2*P3Q2-4.*P1P3*P1Q2*UME2+4.*P1P3*P3Q2*
     . P3Q1-4.*P1P3*P3Q2*EME2-2.*P1P3*P3Q1*EME2+2.*S2P1**2*P3Q2+
     . 2.*S2P1**2*P3Q1+4.*S2P1*S2P3*P1Q2+6.*S2P1*S2P3*P3Q2-4.*
     . S2P1*S2P3*P3Q1+4.*S2P1*S2P3*EME2+2.*S2P1*S2Q2*P3Q2+2.*
     . S2P1*S2Q2*P3Q1+2.*S2P1*S2Q1*P3Q2+2.*S2P1*S2Q1*P3Q1-8.*
     . S2P3**2*P1Q2+4.*S2P3*S2Q2*P1Q2-2.*S2P3*S2Q1*P1Q2+2.*S2P3*
     . S2Q1*EME2-4.*S22*P1Q2*P3Q2-S22*P1Q2*P3Q1+6.*P1Q2**2*P3Q2+
     . 2.*P1Q2**2*P3Q1+6.*P1Q2**2*UME2-2.*P1Q2*P3Q2*P3Q1+P1Q2*
     . P3Q1**2-P1Q2*P3Q1*EME2-P1Q2*P3Q1*UME2+2.*P3Q2**2*EME2)+
     . ANS12+ANS15
      ANS10=4.*Q22**2*(4.*P1P3**3-4.*P1P3**2*P1Q2+4.*P1P3**2*
     . P3Q2+8.*P1P3**2*P3Q1+4.*P1P3*S2P1*S2P3-2.*P1P3*S2P1*S2Q2-
     . 2.*P1P3*S2P1*S2Q1+4.*P1P3*S2P3**2-2.*P1P3*S2P3*S2Q2-2.*
     . P1P3*S22*P1Q2+P1P3*S22*P3Q2+2.*P1P3*P1Q2**2-8.*P1P3*P1Q2*
     . P3Q2-8.*P1P3*P1Q2*P3Q1-4.*P1P3*P1Q2*UME2+4.*P1P3*P3Q2*
     . P3Q1+4.*P1P3*P3Q2*EME2+4.*P1P3*P3Q1**2+2.*P1P3*P3Q1*EME2+
     . 2.*S2P1**2*P3Q2+2.*S2P1**2*P3Q1+4.*S2P1*S2P3*P1Q2+2.*S2P1*
     . S2P3*P3Q2+2.*S2P1*S2P3*P3Q1-2.*S2P1*S2P3*EME2-2.*S2P1*
     . S2P3*UME2+S2P1*S2Q2*UME2-4.*S2P3**2*P1Q2+2.*S2P3*S2Q2*
     . P1Q2-S2P3*S2Q2*EME2+2.*S2P3*S2Q1*P1Q2-2.*S2P3*S2Q1*EME2-
     . 2.*S22*P1Q2*P3Q2-2.*S22*P1Q2*P3Q1+6.*P1Q2**2*P3Q2+4.*P1Q2
     . **2*P3Q1+6.*P1Q2**2*UME2-2.*P1Q2*P3Q2*P3Q1-3.*P1Q2*P3Q1**
     . 2-P1Q2*P3Q1*EME2-P1Q2*P3Q1*UME2+2.*P3Q2**2*EME2+4.*P3Q2*
     . P3Q1*EME2+2.*P3Q1**2*EME2)+2.*Q22*Q32**2*(2.*P1P3**2-2.*
     . P1P3*P1Q2-2.*P1P3*P3Q2-P1P3*P3Q1+P1P3*EME2-P1P3*UME2+2.*
     . S2P1*S2P3+2.*P1Q2*P3Q2+2.*P1Q2*UME2)+ANS11
      ANS9=2.*Q22**3*(-2.*P1P3**2+2.*P1P3*P1Q2-2.*P1P3*P3Q2-P1P3
     . *P3Q1+P1P3*EME2-P1P3*UME2+2.*P1Q2*P3Q2+2.*P1Q2*P3Q1+2.*
     . P1Q2*UME2)+2.*Q22**2*Q32*(-P1P3*S22-4.*P1P3*P3Q2-2.*P1P3*
     . P3Q1-2.*P1P3*EME2-2.*P1P3*UME2+4.*S2P1*S2P3+4.*P1Q2*P3Q2+
     . 2.*P1Q2*P3Q1+4.*P1Q2*UME2-2.*P3Q1*EME2)+ANS10
      ANS5=2.*Q12*Q32**2*(3.*P1P3*S22-4.*S2P1*S2P3-S22*P3Q2)+4.*
     . Q12*Q32*(-8.*P1P3*S2P1*S2P3+8.*P1P3*S2P3*S2Q2-4.*P1P3*
     . S2Q2**2-4.*P1P3*S2Q2*S2Q1+4.*P1P3*S22*P1Q2-2.*P1P3*S22*
     . P3Q2+6.*S2P1*S2P3*P3Q2-4.*S2P1*S2P3*EME2+4.*S2P1*S2Q2*
     . P3Q2+2.*S2P1*S2Q2*P3Q1+2.*S2P1*S2Q1*P3Q2-8.*S2P3**2*P1Q2+
     . 4.*S2P3*S2Q2*P1Q2-2.*S2P3*S2Q2*P3Q2+2.*S2P3*S2Q2*P3Q1+2.*
     . S2P3*S2Q2*EME2+2.*S2P3*S2Q2*UME2+2.*S2P3*S2Q1*P1Q2-2.*
     . S2Q2**2*P3Q1-S2Q2**2*UME2+2.*S2Q2*S2Q1*P3Q2-6.*S22*P1Q2*
     . P3Q2-2.*S22*P1Q2*P3Q1+S22*P3Q2**2)+ANS6+ANS9
      ANS1=2.*Q12*Q22**2*(-P1P3*S22-4.*P1P3*P3Q2-4.*P1P3*P3Q1-2.
     . *P1P3*EME2-4.*P1P3*UME2+4.*S2P1*S2P3-4.*S2P3**2+2.*S2P3*
     . S2Q2+2.*S2P3*S2Q1-2.*S22*P3Q2-2.*S22*P3Q1+8.*P1Q2*P3Q2+5.
     . *P1Q2*P3Q1+8.*P1Q2*UME2-P3Q2*P3Q1-P3Q2*EME2-P3Q2*UME2-2.*
     . P3Q1**2-2.*P3Q1*EME2-2.*P3Q1*UME2)-Q12*Q22*Q32**2*P1P3+2.
     . *Q12*Q22*Q32*(4.*P1P3**2-4.*P1P3*P1Q2-4.*P1P3*P3Q2+2.*
     . P1P3*EME2+2.*S2P1*S2P3-4.*S2P3**2+2.*S2P3*S2Q2-3.*S22*
     . P3Q2-2.*S22*P3Q1+4.*P1Q2*P3Q2+P1Q2*P3Q1+4.*P1Q2*UME2-P3Q2
     . *P3Q1-P3Q2*EME2-P3Q2*UME2)+ANS2+ANS5
      TAU4=Q12**3*Q22*(-P1P3+P3Q2)+2.*Q12**3*S22*(P1P3-P3Q2)+Q12
     . **2*Q22**2*(3.*P3Q2+2.*P3Q1+2.*UME2)+Q12**2*Q22*Q32*(-2.*
     . P1P3+P3Q2)+2.*Q12**2*Q22*(2.*P1P3**2-2.*P1P3*P1Q2-2.*P1P3
     . *P3Q2+P1P3*P3Q1+P1P3*EME2+P1P3*UME2-4.*S2P3**2+2.*S2P3*
     . S2Q2-3.*S22*P3Q2-2.*S22*P3Q1+4.*P1Q2*P3Q2+P1Q2*P3Q1+2.*
     . P1Q2*UME2-P3Q2*P3Q1-P3Q2*EME2-P3Q2*UME2)+2.*Q12**2*Q32*(
     . 3.*P1P3*S22-2.*S2P1*S2P3-2.*S22*P3Q2)+4.*Q12**2*(-4.*P1P3*
     . S2P1*S2P3+4.*P1P3*S2P3*S2Q2-2.*P1P3*S2Q2**2-2.*P1P3*S2Q2*
     . S2Q1+2.*P1P3*S22*P1Q2-P1P3*S22*P3Q2+2.*S2P1*S2P3*P3Q2-2.*
     . S2P1*S2P3*P3Q1-2.*S2P1*S2P3*EME2-2.*S2P1*S2P3*UME2+2.*
     . S2P1*S2Q2*P3Q2+2.*S2P1*S2Q2*P3Q1+S2P1*S2Q2*UME2-4.*S2P3**
     . 2*P1Q2+2.*S2P3*S2Q2*P1Q2-2.*S2P3*S2Q2*P3Q2+2.*S2P3*S2Q2*
     . P3Q1+S2P3*S2Q2*EME2+2.*S2P3*S2Q2*UME2-2.*S2Q2**2*P3Q1-
     . S2Q2**2*UME2+2.*S2Q2*S2Q1*P3Q2-4.*S22*P1Q2*P3Q2-S22*P1Q2*
     . P3Q1+S22*P3Q2**2)+Q12*Q22**3*(P1P3+2.*P3Q2+2.*P3Q1+2.*
     . UME2)+2.*Q12*Q22**2*Q32*(P3Q2+P3Q1+UME2)+ANS1
      ANS5=-P2Q1*P3Q2*EME2+3.*P2Q1*P3Q2*UME2-2.*P2Q1*P3Q1**2-2.*
     . P2Q1*P3Q1*EME2+2.*P2Q1*P3Q1*UME2+2.*P3Q2**2*EME2
      ANS4=4.*P2P3**3-4.*P2P3**2*P2Q1+4.*P2P3**2*P3Q2-4.*P2P3**2
     . *P3Q1+8.*P2P3*S1P2**2+10.*P2P3*S1P2*S1Q2+2.*P2P3*S1P2*
     . S1Q1-4.*P2P3*S1P3**2+4.*P2P3*S1P3*S1Q2+2.*P2P3*S1P3*S1Q1+
     . 4.*P2P3*S1Q2**2+2.*P2P3*S1Q2*S1Q1-2.*P2P3*S12*P3Q2-2.*P2P3
     . *S12*P3Q1+2.*P2P3*P2Q1**2-4.*P2P3*P2Q1*P3Q2+4.*P2P3*P2Q1*
     . P3Q1+4.*P2P3*P2Q1*UME2-4.*P2P3*P3Q2*P3Q1+4.*P2P3*P3Q2*
     . EME2+2.*P2P3*P3Q1*EME2+6.*S1P2**2*P3Q2+2.*S1P2**2*P3Q1-4.
     . *S1P2*S1P3*P2Q1-6.*S1P2*S1P3*P3Q2-2.*S1P2*S1P3*P3Q1-4.*
     . S1P2*S1P3*EME2-4.*S1P2*S1P3*UME2+8.*S1P2*S1Q2*P3Q2+4.*
     . S1P2*S1Q2*P3Q1+4.*S1P2*S1Q2*UME2+4.*S1P2*S1Q1*P3Q2+2.*
     . S1P2*S1Q1*P3Q1+2.*S1P2*S1Q1*UME2-2.*S1P3*S1Q2*P2Q1-4.*
     . S1P3*S1Q2*P3Q2-2.*S1P3*S1Q2*P3Q1-2.*S1P3*S1Q2*EME2-4.*
     . S1P3*S1Q2*UME2-2.*S1P3*S1Q1*P3Q2-2.*S1P3*S1Q1*P3Q1-2.*
     . S1P3*S1Q1*UME2+4.*S1Q2**2*P3Q2+2.*S1Q2**2*P3Q1+4.*S1Q2**2
     . *UME2+4.*S1Q2*S1Q1*P3Q2+2.*S1Q2*S1Q1*P3Q1+4.*S1Q2*S1Q1*
     . UME2+S1Q1**2*UME2-5.*S12*P2Q1*P3Q2-2.*S12*P2Q1*P3Q1+S12*
     . P3Q1**2+4.*P2Q1**2*P3Q2+3.*P2Q1*P3Q2*P3Q1+ANS5
      ANS3=4.*Q12*Q22*ANS4
      ANS9=4.*S1P2*S1Q2*P3Q2*P3Q1-4.*S1P2*S1Q2*P3Q2*EME2-2.*S1P2
     . *S1Q2*P3Q1*EME2+2.*S1P2*S1Q1*P2Q1*P3Q2+2.*S1P2*S1Q1*P3Q2
     . **2+2.*S1P2*S1Q1*P3Q2*P3Q1-2.*S1P2*S1Q1*P3Q2*EME2-2.*S1P3
     . *S1Q2*P2Q1*P3Q2+2.*S1P3*S1Q2*P2Q1*EME2-2.*S1P3*S1Q2*P2Q1*
     . UME2+2.*S1P3*S1Q1*P2Q1*EME2+2.*S1Q2**2*P2Q1*P3Q2+2.*S1Q2
     . **2*P2Q1*UME2-4.*S1Q2**2*P3Q2*EME2-2.*S1Q2**2*P3Q1*EME2+
     . 2.*S1Q2*S1Q1*P2Q1*P3Q2+S1Q2*S1Q1*P2Q1*UME2-4.*S1Q2*S1Q1*
     . P3Q2*EME2-S1Q2*S1Q1*P3Q1*EME2-S1Q1**2*P3Q2*EME2-2.*S12*
     . P2Q1**2*P3Q2-S12*P2Q1*P3Q2*P3Q1+2.*P2Q1**3*P3Q2+2.*P2Q1**
     . 2*P3Q2*P3Q1+2.*P2Q1**2*P3Q2*UME2-2.*P2Q1*P3Q2*P3Q1*EME2
      ANS8=8.*P2P3**2*S1P2*S1P3-8.*P2P3**2*S1P2*S1Q2-4.*P2P3**2*
     . S1P2*S1Q1+8.*P2P3**2*S1P3*S1Q2+4.*P2P3**2*S1P3*S1Q1-8.*
     . P2P3**2*S1Q2**2-8.*P2P3**2*S1Q2*S1Q1-2.*P2P3**2*S1Q1**2-
     . 4.*P2P3**2*P2Q1*P3Q2-8.*P2P3**2*P2Q1*P3Q1+8.*P2P3*S1P2**2*
     . P3Q2+4.*P2P3*S1P2**2*P3Q1+8.*P2P3*S1P2*S1P3*P3Q2+4.*P2P3*
     . S1P2*S1P3*P3Q1+4.*P2P3*S1P2*S1Q2*P3Q2+2.*P2P3*S1P2*S1Q2*
     . P3Q1-4.*P2P3*S1P2*S1Q1*P2Q1+2.*P2P3*S1P2*S1Q1*P3Q2+2.*
     . P2P3*S1P2*S1Q1*P3Q1-8.*P2P3*S1P3**2*P2Q1+8.*P2P3*S1P3*
     . S1Q2*P2Q1+8.*P2P3*S1P3*S1Q2*P3Q2+4.*P2P3*S1P3*S1Q2*P3Q1+
     . 4.*P2P3*S1P3*S1Q1*P2Q1+4.*P2P3*S1P3*S1Q1*P3Q2-4.*P2P3*S1Q2
     . **2*P3Q2-2.*P2P3*S1Q2**2*P3Q1-2.*P2P3*S1Q2*S1Q1*P2Q1-4.*
     . P2P3*S1Q2*S1Q1*P3Q2-2.*P2P3*S1Q1**2*P3Q2+2.*P2P3*S12*P2Q1
     . **2-4.*P2P3*S12*P2Q1*P3Q2-4.*P2P3*S12*P2Q1*P3Q1+4.*P2P3*
     . P2Q1**2*P3Q1+4.*P2P3*P2Q1**2*UME2-4.*P2P3*P2Q1*P3Q2**2-8.
     . *P2P3*P2Q1*P3Q2*P3Q1+4.*S1P2**2*P2Q1*P3Q2+4.*S1P2**2*P3Q2
     . **2+4.*S1P2**2*P3Q2*P3Q1-4.*S1P2*S1P3*P2Q1*P3Q2+4.*S1P2*
     . S1Q2*P2Q1*P3Q2+4.*S1P2*S1Q2*P3Q2**2+ANS9
      ANS7=8.*Q12*ANS8
      ANS13=-2.*S1P3*S1Q2*P2Q1*P3Q1-2.*S1P3*S1Q2*P2Q1*EME2-2.*
     . S1P3*S1Q2*P2Q1*UME2-2.*S1P3*S1Q1*P2Q1*P3Q2-4.*S1P3*S1Q1*
     . P2Q1*P3Q1-2.*S1P3*S1Q1*P2Q1*EME2-4.*S1P3*S1Q1*P2Q1*UME2+
     . 2.*S1Q2**2*P2Q1*P3Q2+2.*S1Q2**2*P2Q1*P3Q1+2.*S1Q2**2*P2Q1*
     . UME2+6.*S1Q2*S1Q1*P2Q1*P3Q2+2.*S1Q2*S1Q1*P2Q1*P3Q1+5.*
     . S1Q2*S1Q1*P2Q1*UME2-2.*S1Q2*S1Q1*P3Q2*EME2-S1Q2*S1Q1*P3Q1
     . *EME2+2.*S1Q1**2*P2Q1*UME2-S1Q1**2*P3Q2*EME2-4.*S12*P2Q1
     . **2*P3Q2-2.*S12*P2Q1*P3Q2**2-S12*P2Q1*P3Q2*P3Q1+2.*S12*
     . P2Q1*P3Q1**2
      ANS12=-8.*P2P3**2*S1P2*S1P3+8.*P2P3**2*S1P2*S1Q2+4.*P2P3**
     . 2*S1P2*S1Q1+4.*P2P3**2*S1P3*S1Q1-4.*P2P3**2*S1Q2*S1Q1-2.*
     . P2P3**2*S1Q1**2-8.*P2P3*S1P2**2*P3Q2-4.*P2P3*S1P2**2*P3Q1
     . -8.*P2P3*S1P2*S1P3*P3Q2+4.*P2P3*S1P2*S1P3*P3Q1+4.*P2P3*
     . S1P2*S1Q2*P2Q1+4.*P2P3*S1P2*S1Q2*P3Q2-6.*P2P3*S1P2*S1Q2*
     . P3Q1+8.*P2P3*S1P2*S1Q1*P2Q1+10.*P2P3*S1P2*S1Q1*P3Q2+2.*
     . P2P3*S1P2*S1Q1*P3Q1-8.*P2P3*S1P3**2*P2Q1+8.*P2P3*S1P3*
     . S1Q2*P2Q1+4.*P2P3*S1P3*S1Q1*P2Q1+4.*P2P3*S1P3*S1Q1*P3Q2+
     . 2.*P2P3*S1Q2**2*P2Q1+4.*P2P3*S1Q2*S1Q1*P2Q1-2.*P2P3*S1Q2*
     . S1Q1*P3Q2+2.*P2P3*S1Q2*S1Q1*P3Q1-2.*P2P3*S1Q1**2*P3Q2-2.*
     . P2P3*S12*P2Q1**2-4.*P2P3*S12*P2Q1*P3Q2-4.*P2P3*S12*P2Q1*
     . P3Q1+4.*S1P2**2*P2Q1*P3Q1-4.*S1P2**2*P3Q2**2-4.*S1P2*S1P3
     . *P2Q1**2-4.*S1P2*S1P3*P2Q1*P3Q2+4.*S1P2*S1P3*P3Q2*P3Q1+2.
     . *S1P2*S1Q2*P2Q1*P3Q2+4.*S1P2*S1Q2*P2Q1*P3Q1-2.*S1P2*S1Q2*
     . P3Q2*P3Q1+4.*S1P2*S1Q2*P3Q2*EME2+2.*S1P2*S1Q2*P3Q1*EME2+
     . 8.*S1P2*S1Q1*P2Q1*P3Q2+4.*S1P2*S1Q1*P3Q2**2+2.*S1P2*S1Q1*
     . P3Q2*EME2-4.*S1P3**2*P2Q1*P3Q2-2.*S1P3*S1Q2*P2Q1**2+ANS13
      ANS11=8.*Q22*ANS12
      ANS17=2.*S1P3*S1Q1*P2Q1*EME2-2.*S1Q2**2*P2Q1*P3Q2-2.*S1Q2
     . **2*P2Q1*UME2-S1Q2*S1Q1*P2Q1*UME2-2.*S1Q2*S1Q1*P3Q2*EME2-
     . S1Q2*S1Q1*P3Q1*EME2-S1Q1**2*P3Q2*EME2-2.*S12*P2Q1**2*P3Q2
     . -2.*S12*P2Q1*P3Q2**2-S12*P2Q1*P3Q2*P3Q1
      ANS16=-8.*P2P3**2*S1P2*S1P3+8.*P2P3**2*S1P2*S1Q2+4.*P2P3**
     . 2*S1P2*S1Q1+4.*P2P3**2*S1P3*S1Q1-4.*P2P3**2*S1Q2*S1Q1-2.*
     . P2P3**2*S1Q1**2-8.*P2P3*S1P2**2*P3Q2-4.*P2P3*S1P2**2*P3Q1
     . -8.*P2P3*S1P2*S1P3*P3Q2+4.*P2P3*S1P2*S1P3*P3Q1-4.*P2P3*
     . S1P2*S1Q2*P2Q1+4.*P2P3*S1P2*S1Q2*P3Q2-6.*P2P3*S1P2*S1Q2*
     . P3Q1-4.*P2P3*S1P2*S1Q1*P2Q1+10.*P2P3*S1P2*S1Q1*P3Q2+2.*
     . P2P3*S1P2*S1Q1*P3Q1-8.*P2P3*S1P3**2*P2Q1+8.*P2P3*S1P3*
     . S1Q2*P2Q1+4.*P2P3*S1P3*S1Q1*P2Q1+4.*P2P3*S1P3*S1Q1*P3Q2-
     . 2.*P2P3*S1Q2**2*P2Q1-2.*P2P3*S1Q2*S1Q1*P2Q1-2.*P2P3*S1Q2*
     . S1Q1*P3Q2+2.*P2P3*S1Q2*S1Q1*P3Q1-2.*P2P3*S1Q1**2*P3Q2+2.*
     . P2P3*S12*P2Q1**2-4.*P2P3*S12*P2Q1*P3Q2-4.*P2P3*S12*P2Q1*
     . P3Q1-4.*S1P2**2*P3Q2**2-4.*S1P2*S1P3*P2Q1*P3Q2+4.*S1P2*
     . S1P3*P3Q2*P3Q1-2.*S1P2*S1Q2*P2Q1*P3Q2-2.*S1P2*S1Q2*P3Q2*
     . P3Q1+4.*S1P2*S1Q2*P3Q2*EME2+2.*S1P2*S1Q2*P3Q1*EME2+2.*
     . S1P2*S1Q1*P2Q1*P3Q2+4.*S1P2*S1Q1*P3Q2**2+2.*S1P2*S1Q1*
     . P3Q2*EME2-4.*S1P3**2*P2Q1*P3Q2+4.*S1P3*S1Q2*P2Q1*P3Q2+2.*
     . S1P3*S1Q2*P2Q1*EME2+2.*S1P3*S1Q2*P2Q1*UME2+2.*S1P3*S1Q1*
     . P2Q1*P3Q2+ANS17
      ANS15=8.*Q32*ANS16
      ANS18=32.*P2Q1*(4.*P2P3**2*S1P3*S1Q2+4.*P2P3**2*S1P3*S1Q1-
     . 4.*P2P3**2*S1Q2**2-6.*P2P3**2*S1Q2*S1Q1-2.*P2P3**2*S1Q1**2
     . +4.*P2P3*S1P2*S1P3*P3Q1+4.*P2P3*S1P2*S1Q2*P3Q2-2.*P2P3*
     . S1P2*S1Q2*P3Q1+6.*P2P3*S1P2*S1Q1*P3Q2+2.*P2P3*S1P2*S1Q1*
     . P3Q1-4.*P2P3*S1P3**2*P2Q1+4.*P2P3*S1P3*S1Q2*P2Q1+4.*P2P3*
     . S1P3*S1Q2*P3Q2+2.*P2P3*S1P3*S1Q2*P3Q1+2.*P2P3*S1P3*S1Q1*
     . P2Q1+4.*P2P3*S1P3*S1Q1*P3Q2-2.*P2P3*S1Q2**2*P3Q2-P2P3*
     . S1Q2**2*P3Q1-3.*P2P3*S1Q2*S1Q1*P3Q2+P2P3*S1Q2*S1Q1*P3Q1-
     . 2.*P2P3*S1Q1**2*P3Q2-2.*P2P3*S12*P2Q1*P3Q2-2.*P2P3*S12*
     . P2Q1*P3Q1+2.*S1P2**2*P3Q2*P3Q1-2.*S1P2*S1P3*P2Q1*P3Q2+2.*
     . S1P2*S1P3*P3Q2*P3Q1+2.*S1P2*S1Q2*P3Q2**2+S1P2*S1Q2*P3Q2*
     . P3Q1+2.*S1P2*S1Q1*P2Q1*P3Q2+3.*S1P2*S1Q1*P3Q2**2+S1P2*
     . S1Q1*P3Q2*P3Q1-2.*S1P3**2*P2Q1*P3Q2+S1P3*S1Q2*P2Q1*P3Q2+
     . S1P3*S1Q1*P2Q1*P3Q2-2.*S1Q2**2*P3Q2*EME2-S1Q2**2*P3Q1*
     . EME2+S1Q2*S1Q1*P2Q1*P3Q2-3.*S1Q2*S1Q1*P3Q2*EME2-S1Q2*S1Q1
     . *P3Q1*EME2-S1Q1**2*P3Q2*EME2-S12*P2Q1**2*P3Q2-S12*P2Q1*
     . P3Q2**2-S12*P2Q1*P3Q2*P3Q1)
      ANS14=4.*Q32**2*(4.*P2P3*S1P2**2+2.*P2P3*S1P2*S1Q2-2.*P2P3
     . *S1P2*S1Q1+2.*P2P3*S12*P2Q1-P2P3*S12*P3Q1+2.*S1P2**2*P3Q2
     . -2.*S1P2*S1P3*P3Q2-2.*S1P2*S1P3*EME2-2.*S1P2*S1P3*UME2+2.
     . *S1P2*S1Q2*P3Q2+2.*S1P2*S1Q2*UME2+S1P2*S1Q1*UME2+S1P3*
     . S1Q1*EME2)+ANS15+ANS18
      ANS10=4.*Q22**2*(-4.*P2P3*S1P2**2-2.*P2P3*S1P2*S1Q2+4.*
     . P2P3*S1P2*S1Q1+2.*P2P3*S1Q2*S1Q1-2.*P2P3*S12*P2Q1-P2P3*
     . S12*P3Q1-2.*S1P2**2*P3Q2+2.*S1P2**2*P3Q1-4.*S1P2*S1P3*
     . P2Q1+2.*S1P2*S1P3*P3Q2+2.*S1P2*S1P3*P3Q1+2.*S1P2*S1P3*
     . EME2+2.*S1P2*S1P3*UME2-2.*S1P2*S1Q2*P3Q2-2.*S1P2*S1Q2*
     . UME2+2.*S1P2*S1Q1*P3Q2-2.*S1P2*S1Q1*P3Q1-S1P2*S1Q1*UME2-
     . 2.*S1P3*S1Q2*P2Q1-2.*S1P3*S1Q1*P3Q2-2.*S1P3*S1Q1*P3Q1-S1P3
     . *S1Q1*EME2-2.*S1P3*S1Q1*UME2+2.*S1Q2*S1Q1*P3Q2+2.*S1Q2*
     . S1Q1*UME2+S1Q1**2*UME2-2.*S12*P2Q1*P3Q2+S12*P3Q1**2)+4.*
     . Q22*Q32*(2.*P2P3*S1P2*S1Q1+2.*P2P3*S1Q2*S1Q1-2.*P2P3*S12*
     . P3Q1+2.*S1P2**2*P3Q1-4.*S1P2*S1P3*P2Q1+2.*S1P2*S1P3*P3Q1+
     . 2.*S1P2*S1Q1*P3Q2-2.*S1P2*S1Q1*P3Q1-2.*S1P3*S1Q2*P2Q1-2.*
     . S1P3*S1Q1*P3Q2-2.*S1P3*S1Q1*P3Q1-2.*S1P3*S1Q1*UME2+2.*
     . S1Q2*S1Q1*P3Q2+2.*S1Q2*S1Q1*UME2+S1Q1**2*UME2-2.*S12*P2Q1
     . *P3Q2+S12*P3Q1**2)+ANS11+ANS14
      ANS6=2.*Q12*Q32**2*P2P3*(2.*P2P3+S12-2.*P2Q1+2.*P3Q2+P3Q1-
     . EME2+UME2)+4.*Q12*Q32*(4.*P2P3**3+4.*P2P3**2*P2Q1+4.*P2P3
     . **2*P3Q2-4.*P2P3**2*P3Q1-2.*P2P3*S1P2*S1Q2-4.*P2P3*S1P2*
     . S1Q1-4.*P2P3*S1P3**2+4.*P2P3*S1P3*S1Q2+2.*P2P3*S1P3*S1Q1-
     . 2.*P2P3*S1Q2*S1Q1+4.*P2P3*S12*P2Q1-2.*P2P3*S12*P3Q2-2.*
     . P2P3*S12*P3Q1-2.*P2P3*P2Q1**2+4.*P2P3*P2Q1*P3Q2+4.*P2P3*
     . P2Q1*P3Q1+4.*P2P3*P2Q1*UME2-4.*P2P3*P3Q2*P3Q1+4.*P2P3*
     . P3Q2*EME2+2.*P2P3*P3Q1*EME2+2.*S1P2**2*P3Q2-2.*S1P2*S1P3*
     . P3Q2+2.*S1P3*S1Q2*EME2+2.*S1P3*S1Q1*EME2-S12*P2Q1*P3Q2+2.
     . *P2Q1**2*P3Q2+P2Q1*P3Q2*P3Q1+P2Q1*P3Q2*EME2+P2Q1*P3Q2*
     . UME2+2.*P3Q2**2*EME2)+ANS7+ANS10
      ANS2=2.*Q12*Q22*Q32*(4.*P2P3*P3Q1-4.*S1P2*S1P3-2.*S1P3*
     . S1Q2-2.*S12*P3Q2-S12*P3Q1-4.*P2Q1*UME2+2.*P3Q2*P3Q1-P3Q1
     . **2-P3Q1*EME2+P3Q1*UME2)+ANS3+ANS6
      ANS1=4.*Q12**2*(-4.*P2P3**3+4.*P2P3**2*P2Q1-8.*P2P3**2*
     . P3Q2-4.*P2P3**2*P3Q1-4.*P2P3*S1P2**2-4.*P2P3*S1P2*S1Q2-2.
     . *P2P3*S1P2*S1Q1-4.*P2P3*S1P3**2+4.*P2P3*S1P3*S1Q2+2.*P2P3
     . *S1P3*S1Q1-2.*P2P3*S1Q2**2-2.*P2P3*S1Q2*S1Q1+2.*P2P3*S12*
     . P2Q1-2.*P2P3*S12*P3Q2-P2P3*S12*P3Q1-2.*P2P3*P2Q1**2+4.*
     . P2P3*P2Q1*P3Q2+4.*P2P3*P2Q1*P3Q1+4.*P2P3*P2Q1*UME2-4.*
     . P2P3*P3Q2**2-4.*P2P3*P3Q2*P3Q1-4.*P2P3*P3Q2*EME2-2.*P2P3*
     . P3Q1*EME2+2.*S1P2*S1P3*EME2+2.*S1P2*S1P3*UME2-2.*S1P2*
     . S1Q2*UME2-S1P2*S1Q1*UME2+2.*S1P3*S1Q2*EME2+S1P3*S1Q1*EME2
     . -S12*P2Q1*P3Q2+2.*P2Q1**2*P3Q2+P2Q1*P3Q2*P3Q1-P2Q1*P3Q2*
     . EME2+P2Q1*P3Q2*UME2-2.*P3Q2**2*EME2-2.*P3Q2*P3Q1*EME2)-(
     . 2.*Q12*Q22**2*Q32)*(P2P3+P3Q2+UME2)+2.*Q12*Q22**2*(-2.*
     . P2P3**2-P2P3*S12+2.*P2P3*P2Q1-2.*P2P3*P3Q2+3.*P2P3*P3Q1+
     . P2P3*EME2-P2P3*UME2-4.*S1P2*S1P3-2.*S1P3*S1Q2-2.*S12*P3Q2
     . -S12*P3Q1+2.*P2Q1*P3Q2+2.*P3Q2*P3Q1-P3Q1**2-P3Q1*EME2+
     . P3Q1*UME2)-(2.*Q12*Q22*Q32**2)*(P2P3+P3Q2+UME2)+ANS2
      TAU5=-Q12**4*P2P3+Q12**3*Q22*(2.*P3Q2+P3Q1)-2.*Q12**3*Q32*
     . P2P3+2.*Q12**3*(2.*P2P3**2+P2P3*S12-2.*P2P3*P2Q1+2.*P2P3*
     . P3Q2+P2P3*P3Q1+P2P3*EME2+P2P3*UME2+P2Q1*P3Q2)+Q12**2*Q22
     . **2*(P2P3+2.*P3Q2+P3Q1)+Q12**2*Q22*Q32*(2.*P2P3+2.*P3Q2+
     . P3Q1-2.*UME2)+2.*Q12**2*Q22*(-2.*P2P3*EME2+4.*P2P3*UME2-
     . 2.*S12*P3Q2-S12*P3Q1+5.*P2Q1*P3Q2+2.*P2Q1*P3Q1-2.*P3Q2*
     . EME2+2.*P3Q2*UME2-P3Q1**2-P3Q1*EME2+P3Q1*UME2)-Q12**2*Q32
     . **2*P2P3+2.*Q12**2*Q32*(4.*P2P3**2+2.*P2P3*S12-4.*P2P3*
     . P2Q1+4.*P2P3*P3Q2+2.*P2P3*P3Q1+2.*P2P3*UME2+P2Q1*P3Q2)+
     . ANS1
 
 
      RETURN
      END

