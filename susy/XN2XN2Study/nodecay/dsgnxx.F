C* (Update Record)
C*    95/06/24  K.Fujii          Original version.
C*
 
      SUBROUTINE DSGNXX(IX1,IX2,RS,IHEL,CSX,FIX,SG)
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SSCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
      PARAMETER    ( IxPRC = 2 )
      INTEGER*4    IX1, IX2, IHEL(4)
      REAL   *8    RS, CSX, FIX, SG
C--
      REAL   *8    PV(0:3,4), QV(0:3,4), AMP2(0:2), EC(3,3), EB(3,3)
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C  Spin average is taken only for positron.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         NB     = 1
         FACT   = xGV2FB*x2PI/(x4PI**(3*NB))
         AME    = AMF(1,2,1)
         AM1    = SZM(IX1)
         AM2    = SZM(IX2)
         CALL UVZERO(9*IxPRC,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(16*IxPRC,PV)
      CALL UVZERO(16*IxPRC,QV)
C
C-- e+ + e- --> CM.
C
      S       = RS*RS
      EBM     = RS/2
      PBM     = SQRT((EBM-AME)*(EBM+AME))
      BT0     = PBM/EBM
      PV(0,1) = EBM
      PV(3,1) = PBM
      CALL PMIRRx(PV(0,1),PV(0,2))
      QV(0,1) = RS
C
C-- CM --> x x-bar
C
      AM12 = AM1*AM1
      AM22 = AM2*AM2
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,CSX,FIX,
     .            PV(0,3),PV(0,4),BT1,EB)
      IF ( BT1.EQ.0. )                           GO TO 9999
C--
C  Calculate phase space weight.
C--
      IF (IX1 .EQ. IX2) THEN
         IDNT = 2
      ELSE
         IDNT = 1
      ENDIF
      WAT = FACT*BT1/(2*S*BT0)/IDNT
C--
C  Calculate amplitude squared.
C--
      CALL FULNXX(IX1,IX2,IHEL,PV,AMP2)
C--
C  Differenctial cross section.
C--
      SG  = AMP2(0)*WAT
C--
C  That's it.
C--
      RETURN
C--
C  Kinematically forbidden.
C--
9999  SG = 0
      RETURN
      END
