C************************************************
CC********************************************************************CC
C*                                                                    *C
C*===========-----=======                                             *C
C*  Subroutine U2VANG(A,B,TH,PHI)                                     *C
C*===========-----=======                                             *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    Calculate angle between two vector.
C* (Inputs)                                                           *C
C*      A(0:3), b(0:3); 2 3-vector.                                   *C
C* (Output)                                                           *C
C*      TH, Phi  : Angle in radian
C* (Update Record)                                                    *C
C*  10/27/90  A. Miyamoto  Original version.
C*                                                                    *C
CC********************************************************************CC
C
C =============-----=======
      SUBROUTINE U2VANG(A,B,TH,PHI)
C =============-----=======
C
      DIMENSION A(3),B(3)
 
      PTA = SQRT( A(1)**2 + A(2)**2 )
      PTB = SQRT( B(1)**2 + B(2)**2 )
      TANTHA = PTA/A(3)
      TANTHB = PTB/B(3)
CC    TH = ATAN2( TANTHA-TANTHB, 1. + TANTHA*TANTHB )
C
      TH  = ATAN2( PTA*B(3)-PTB*A(3),  A(3)*B(3) + PTA*PTB )
      PHI = ATAN2( A(2)*B(1)-A(1)*B(2), A(1)*B(1)+A(2)*B(2) )
 
      RETURN
      END
