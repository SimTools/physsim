CC********************************************************************CC
C*                                                                    *C
C*==========================--===                                     *C
C*  Subroutine UBSTBD(PB,PR,PA)                                       *C
C*==========================--===                                     *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    routine to transform PB(4) in PR(4)-rest frame to PA(4)         *C
C*    in PR(4)-moving frame.                                          *C
C* (Inputs)                                                           *C
C*     PB(4) = 4-vector of a particle in PR-rest frame.               *C
C*     PR(4) = 4-vector of a reference particle.                      *C
C* (Output)                                                           *C
C*     PA(4) = 4-vector of a particle in PR-moving frame.             *C
C* (Relation)                                                         *C
C*    Calls                                                           *C
C*       UCNJ4D, UBSTAD                                               *C
C* (Update Record)                                                    *C
C*    7/29/85  K. Fujii  Original version.                            *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UBSTBD(PB,PR,PA)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     PB(4), PR(4), PA(4), PRNV(4)
      CALL UCNJ4D(PR,PRNV)
      CALL UBSTAD(PB,PRNV,PA)
      RETURN
      END
