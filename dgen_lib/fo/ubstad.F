CC********************************************************************CC
C*                                                                    *C
C*==========================--===                                     *C
C*  Subroutine UBSTAD(PB,PR,PA)                                       *C
C*==========================--===                                     *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    routine to transform PB(4) to PR(4)-rest frame ,PA(4).          *C
C* (Inputs)                                                           *C
C*     PB(4) = 4-vector of a particle.                                *C
C*     PR(4) = 4-vector of a reference particle.                      *C
C* (Output)                                                           *C
C*     PA(4) = 4-vector of a particle in PR-rest frame.               *C
C* (Relation)                                                         *C
C*     Calls                                                          *C
C*        UDOT4D                                                      *C
C* (Update Record)                                                    *C
C*    7/29/85  K. Fujii  Original version.                            *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UBSTAD(PB,PR,PA)
 
      IMPLICIT  REAL*8 ( A-H, O-Z )
      REAL*8    PB(4), PR(4), PA(4)
      AM = UDOT4D(PR,PR)
      IF ( AM.LE.0.D0 )                          GO TO 9999
      AM    = SQRT(AM)
      PA4   = UDOT4D(PB,PR)/AM
      A     = (PB(4)+PA4)/(PR(4)+AM)
      PA(4) = PA4
      DO 10 I = 1, 3
         PA(I) = PB(I) - A*PR(I)
10    CONTINUE
      RETURN
9999  WRITE(6,*) ' >>>> Error in UBSTAD. -ve mass squared.>>>'
      WRITE(6,*) '   AM2 = ', AM
      RETURN
      END
