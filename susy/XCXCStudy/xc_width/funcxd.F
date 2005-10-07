C* (Update Record)
C*   95/05/10  K.Fujii	Original version for X+ --> X0 fu + fdbar.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "XCXCCN.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 4, NHEL = 4 )
      INTEGER*4   IPV(3,2)
      REAL   *4   SG, PV(0:3,NP), AM(NP), AMR(3,3)
C--
      REAL   *4   AMSF(3,2,2), GTSF(3,2,2)
C--
      REAL   *8   WAT, WT
C--
      INTEGER*4   IHEL(NP,NHEL,2)
      DATA        (( IHEL(I,J,1), I=1,NP ), J=1,NHEL )  /
     .                   -1,   -1, -1,+1,  
     .                   -1,   +1, -1,+1,  
     .                   -1,   -1, +1,-1,  
     .                   -1,   +1, +1,-1 /
      DATA        (( IHEL(I,J,2), I=1,NP ), J=1,NHEL )  /
     .                   +1,   -1, -1,+1,  
     .                   +1,   +1, -1,+1,  
     .                   +1,   -1, +1,-1,  
     .                   +1,   +1, +1,-1 /
C--
      DATA NCALL /  0 /
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         AMX    = SWM(1)
         GMX    = GMSW(1)
C--
         AMSF(1,1,1) = SFM(1)
         AMSF(1,2,1) = SFM(2)
         AMSF(1,1,2) = SFM(4)
         AMSF(1,2,2) = SFM(6)
C--
         AMSF(2,1,1) = SFM(1)
         AMSF(2,2,1) = SFM(2)
         AMSF(2,1,2) = SFM(4)
         AMSF(2,2,2) = SFM(6)
C--
         AMSF(3,1,1) = SFM(1)
         AMSF(3,2,1) = SFM(2)
         AMSF(3,1,2) = SFM(4)
         AMSF(3,2,2) = SFM(6)
C--
         GTSF(1,1,1) = GMSF(1)
         GTSF(1,2,1) = GMSF(2)
         GTSF(1,1,2) = GMSF(4)
         GTSF(1,2,2) = GMSF(6)
C--
         GTSF(2,1,1) = GMSF(1)
         GTSF(2,2,1) = GMSF(2)
         GTSF(2,1,2) = GMSF(4)
         GTSF(2,2,2) = GMSF(6)
C--
         GTSF(3,1,1) = GMSF(1)
         GTSF(3,2,1) = GMSF(2)
         GTSF(3,1,2) = GMSF(4)
         GTSF(3,2,2) = GMSF(6)
      ENDIF
C--
C  Shuffle integration variables.
C--
      DO 100 I = 1, NZZ
         Z(ISHUFL(I)) = ZZ(I)
100   CONTINUE 
C--
C  Set independent variables.
C     Z( 1) : X+ decay mode
C      ( 2) : helicity combination
C      ( 3) : m(+ab)**2
C      ( 4) : m(+ac)**2
C      ( 5) : cos(theta_a)     in X+ rest frame
C      ( 6) : phi_a            in X+ rest frame
C      ( 7) : phi_b            in X+ rest frame
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  First set initial states.
C--
      IDPTCL(1, 1) = 0
      IDPTCL(2, 1) = 0
      IDPTCL(3, 1) = 0
C--
C  Then select final states.
C--
      IDPTCL(1, 2) = 0
      IDPTCL(2, 2) = 0
      IDPTCL(3, 2) = 0
C-- W+ from X+.
      XW1 = Z(1)
      DO 200 IMD = 1, 12
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1, 3) = MDVDK(1,IMD,1)
      IDPTCL(2, 3) = 1
      IDPTCL(3, 3) = MDVDK(3,IMD,1)
      IDPTCL(1, 4) = MDVDK(2,IMD,1)
      IDPTCL(2, 4) = 2
      IDPTCL(3, 4) = MDVDK(3,IMD,1)
      IMD1         = IMD
C--
      BRW1 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      WAT  = WAT/BRW1
C--
C  Select helicity combination.
C--
      HLM   = 0.5
      IF ( Z(2).LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*Z(2)/HLM + 1
         WAT  = WAT*NHEL
      ELSE
         ICMB = 2
         JCMB = NHEL*(Z(2)-HLM)/(1-HLM) + 1
         WAT  = WAT*NHEL
      ENDIF
      JCMB = MIN(JCMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
C--
C  Select kinematic variables.
C--
      AM(1)   = AMX
      AM(2)   = SZM(1)
      AM(3)   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM(4)   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
C>>>
CTBW       AM(1)   = AMF(3,1,2)
CTBW       AM(2)   = AMF(3,2,2)
C>>>
      IF ( AM(1).LE.AM(2)+AM(3)+AM(4) ) THEN
         FUNC = 0
         RETURN
      ENDIF
C--
      AMR(1,1) = AMSF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3)) 
      AMR(2,1) = GTSF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AMR(3,1) = AMR(1,1) + AM(4)
      AMR(1,2) = AMSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AMR(2,2) = GTSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AMR(3,2) = AMR(1,2) + AM(3)
      AMR(1,3) = AMW
      AMR(2,3) = GMWTOT
      AMR(3,3) = AMR(1,3) + AM(2)
C--
C  Set phase space variables.
C--
      CALL GETXPH(AM,AMR,Z(3),IPV,XPHASE,WT)
      IF ( WT.LE.0.D0 ) THEN
         FUNC = 0
         RETURN
      ENDIF
      WAT = WAT*WT
C--
C  Calculate differential width.
C--
      CALL GMCXDF(IDPTCL,IHLCMB,AM,IPV,XPHASE,SG,PV)
      FUNC   = SG*WAT
      IF ( FUNC.LE.0.D0 ) RETURN
C--
C  Save 4-momenta in PVOUT.
C--
      CALL UVCOPY(4*NP,PV(0,1),PVOUT(0,1))
C--
C  Fill hists. and plots.
C--
      COSLB = PV(3,4)/UABSV3(PV(1,4))
      COSBB = XPHASE(2,1)
      PHIBB = XPHASE(3,1)
      COMB  = JCMB + NHEL*(ICMB-1) + .1
      QXP2  = XPHASE(1,1)
      Q42   = XPHASE(1,2)
      Q52   = XPHASE(2,2)
C--      
      CALL XHFILL(  1, DBLE(COSBB)       ,FUNC )
      CALL XHFILL(  2, DBLE(PHIBB)*xRD2DG,FUNC )
      CALL XHFILL(  3, DBLE(SQRT(Q42))   ,FUNC )
      CALL XHFILL(  4, DBLE(SQRT(Q52))   ,FUNC )
      CALL XHFILL(  5, DBLE(COMB)        ,FUNC )
      IF ( IDPTCL(3,4).EQ.1 .AND. IDPTCL(1,4).LE.2 ) THEN
         CALL XHFILL(  6, DBLE(PV(0,4))     ,FUNC )
         CALL XHFILL(  7, DBLE(COSLB)       ,FUNC )
      ENDIF
      CALL XHFILL(  8, DBLE(IMD1+.1)        ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
