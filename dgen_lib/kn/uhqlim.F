CC**********************************************************************
C*
C*=========================================-----------===
C* Subroutine UHPQLIM(AMPR,AM1,AM2,AM3,Q12,Q13MN,Q13MX)
C*=========================================-----------===
C*
C* (Purpose)
C*    Calculate kinematic boundary for PR ---> 1 + 2 + 3.
C* (Inputs)
C*       AMPR    : (R*16) : parent mass.
C*       AM1     : (R*16) : 1st daughter mass.
C*       AM2     : (R*16) : 2nd daughter mass.
C*       AM3     : (R*16) : 3rd daughter mass.
C*       Q12     : (R*16) : invariant mass for pair 1-2.
C* (Output)
C*       Q13MN   : (R*16) : Q_13 minimum.
C*       Q13MX   : (R*16) : Q_13 maximum.
C* (Relation)
C*    Invokes no external subroutines or functions.
C* (Update Record)
C*    95/04/05  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE UHQLIM(AMPR,AM1,AM2,AM3,Q12,Q13MN,Q13MX)
       
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL*8     AMPR, AM1, AM2, AM3, Q12, Q13MN, Q13MX
C
C========< Entry Point >================================================
C
C--
C  Calculate (E,P) in 1-2 rest frame.
C--
      E1 = Q12*Q12 
      E1 = ( E1 + (AM1-AM2)*(AM1+AM2) )/2/Q12
      P1 = SQRT(MAX(( E1 - AM1 )*( E1 + AM1 ),0.D0))
C--
      E3 = (AMPR-Q12)*(AMPR+Q12)
      E3 = ( E3 - AM3*AM3 )/2/Q12
      P3 = SQRT(MAX(( E3 - AM3 )*( E3 + AM3 ),0.D0))
C--
C  Calculate kinematic boundary.
C--
      E1P3  = E1 + E3
      P1P3  = P1 + P3
      P1M3  = P1 - P3
      Q13MN = SQRT(MAX((E1P3-P1P3)*(E1P3+P1P3),0.D0))
      Q13MX = SQRT(MAX((E1P3-P1M3)*(E1P3+P1M3),0.D0))
C--
C  That's it.
C--
      RETURN
      END
