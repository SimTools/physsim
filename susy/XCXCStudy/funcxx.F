C* (Update Record)
C*   92/06/09  K.Fujii	FUNC for e+ e- --> X+ X- full amplitudes.
C*   95/03/21  K.Fujii	REAL*4 --> REAL*8 for
C*				Q*MN, Q#MX, TH*MN, TH*MX, TH*
C*			with * = 3, 4.
C*  			QED vtx correction ( 1+3*BTE/4 ) and
C*                      ( 1 + 2*ALF0/PI * ( PI**2/6 - 1/4 ) )
C*                      included.
C*   95/04/04  K.Fujii	Added an option to shuffle integration
C*			variables.
C*   95/04/06  K.Fujii	Changed integration variables.
C*   00/06/24  K.Fujii	Modified for Bases++.
C*
      REAL FUNCTION FUNC*8(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "XCXCCN.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 8, NHEL = 16 )
      INTEGER*4   IPV(3,2)
      REAL   *4   SG, PV(0:3,NP), AM(NP), AMX(2), GMX(2), AMR(3,3,2)
C--
      REAL   *4   AMSF(3,2,2), GTSF(3,2,2)
C--
      REAL   *4   QV(0:3,20)
      REAL   *8   WAT, WT
C--
      INTEGER*4   IHEL(NP,NHEL,2)
      DATA        (( IHEL(I,J,1), I=1,NP ), J=1,NHEL )  /
     .                   -1,+1,   -1, -1,+1,   -1, +1,-1,
     .                   -1,+1,   -1, -1,+1,   +1, +1,-1,
     .                   -1,+1,   +1, -1,+1,   -1, +1,-1,
     .                   -1,+1,   +1, -1,+1,   +1, +1,-1,
     .                   -1,+1,   -1, -1,+1,   -1, -1,+1,
     .                   -1,+1,   -1, -1,+1,   +1, -1,+1,
     .                   -1,+1,   +1, -1,+1,   -1, -1,+1,
     .                   -1,+1,   +1, -1,+1,   +1, -1,+1,
     .                   -1,+1,   -1, +1,-1,   -1, +1,-1,
     .                   -1,+1,   -1, +1,-1,   +1, +1,-1,
     .                   -1,+1,   +1, +1,-1,   -1, +1,-1,
     .                   -1,+1,   +1, +1,-1,   +1, +1,-1,
     .                   -1,+1,   -1, +1,-1,   -1, -1,+1,
     .                   -1,+1,   -1, +1,-1,   +1, -1,+1,
     .                   -1,+1,   +1, +1,-1,   -1, -1,+1,
     .                   -1,+1,   +1, +1,-1,   +1, -1,+1/
      DATA        (( IHEL(I,J,2), I=1,NP ), J=1,NHEL )  /
     .                   +1,-1,   -1, -1,+1,   -1, +1,-1,
     .                   +1,-1,   -1, -1,+1,   +1, +1,-1,
     .                   +1,-1,   +1, -1,+1,   -1, +1,-1,
     .                   +1,-1,   +1, -1,+1,   +1, +1,-1,
     .                   +1,-1,   -1, -1,+1,   -1, -1,+1,
     .                   +1,-1,   -1, -1,+1,   +1, -1,+1,
     .                   +1,-1,   +1, -1,+1,   -1, -1,+1,
     .                   +1,-1,   +1, -1,+1,   +1, -1,+1,
     .                   +1,-1,   -1, +1,-1,   -1, +1,-1,
     .                   +1,-1,   -1, +1,-1,   +1, +1,-1,
     .                   +1,-1,   +1, +1,-1,   -1, +1,-1,
     .                   +1,-1,   +1, +1,-1,   +1, +1,-1,
     .                   +1,-1,   -1, +1,-1,   -1, -1,+1,
     .                   +1,-1,   -1, +1,-1,   +1, -1,+1,
     .                   +1,-1,   +1, +1,-1,   -1, -1,+1,
     .                   +1,-1,   +1, +1,-1,   +1, -1,+1/
C--
      DATA NCALL /  0 /
      EXTERNAL   XHFILL
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         AMX(1) = SWM(1)
         GMX(1) = GMSW(1)
         AMX(2) = SWM(1)
         GMX(2) = GMSW(1)
C--
         AME   = AMF(1,2,1)
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
C--
         xALF0 = 1/137.0359895E0
      ENDIF
C--
C  Shuffle integration variables.
C--
      DO 100 I = 1, NZZ
         Z(ISHUFL(I)) = ZZ(I)
100   CONTINUE 
C--
C  Set independent variables.
C     Z( 1) : e- beam
C      ( 2) : e+ beam
C      ( 3) : bremsstrahlng
C      ( 4) : X+ decay mode
C      ( 5) : X- decay mode
C      ( 6) : e- helicity
C      ( 7) : m(+ab)**2
C      ( 8) : m(+ac)**2
C      ( 9) : m(-ab)**2
C      (10) : m(-ac)**2
C      (11) : m(X+)**2
C      (12) : m(X-)**2
C      (13) : cos(theta_X-)
C      (14) : phi_X-
C      (15) : cos(theta_a)     in X+ rest frame
C      (16) : phi_a            in X+ rest frame
C      (17) : phi_b            in X+ rest frame
C      (18) : cos(theta_a)     in X- rest frame
C      (19) : phi_a            in X- rest frame
C      (20) : phi_b            in X- rest frame
C      (21) : E_beam spread of E-
C      (22) : E_beam spread of E+
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Beam energy.
C--
      BMENGY(0) = ROOTS/2
      BMENGY(1) = BMENGY(0)
      BMENGY(2) = BMENGY(0)
      BMENGY(3) = BMENGY(0)
      BMENGY(4) = BMENGY(0)
      DPDEBM    = 1
      RS        = ROOTS
      EMBM      = RS/2
      EPBM      = RS/2
      QED       = 1
C--
C  First decide e- and e+ energy.
C--
      IF ( ISRBM.EQ.3 ) THEN
         EBM = ROOTS/2
         CALL PROBEB(REAL(Z(21)),REAL(Z(1)),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(REAL(Z(22)),REAL(Z(2)),EBM,EPBM,EPBM0,PEP)
         RS  = 2*SQRT(EMBM*EPBM)
         IF ( RS.LT.5. ) THEN
            FUNC = 0
            RETURN
         ENDIF
         BMENGY(0) = EBM
         BMENGY(1) = EMBM0
         BMENGY(2) = EPBM0
         BMENGY(3) = EMBM
         BMENGY(4) = EPBM
         DPDEBM    = PEM*PEP
      ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(1)*1.E4
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
C--
C  First set initial states.
C--
      IDPTCL(1, 1) = 1
      IDPTCL(2, 1) = 2
      IDPTCL(3, 1) = 1
      IDPTCL(1, 2) = 1
      IDPTCL(2, 2) = 2
      IDPTCL(3, 2) = 1
C--
C  Then select final states.
C--
C-- W+ from X+.
      XW1 = Z(4)
      DO 200 IMD = IMD1MN, IMD1MX 
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1, 4) = MDVDK(1,IMD,1)
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = MDVDK(3,IMD,1)
      IDPTCL(1, 5) = MDVDK(2,IMD,1)
      IDPTCL(2, 5) = 2
      IDPTCL(3, 5) = MDVDK(3,IMD,1)
      IMD1         = IMD
C--
      BRW1 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      WAT  = WAT/BRW1
C
C-- W- from X-.
      XW2 = Z(5)
      DO 300 IMD = IMD2MN, IMD2MX
         IF ( XW2.LT.BRVDK(IMD,1) )              GO TO 310
300   CONTINUE
C--
310   IDPTCL(1, 7) = MDVDK(1,IMD,1)
      IDPTCL(2, 7) = 1
      IDPTCL(3, 7) = MDVDK(3,IMD,1)
      IDPTCL(1, 8) = MDVDK(2,IMD,1)
      IDPTCL(2, 8) = 2
      IDPTCL(3, 8) = MDVDK(3,IMD,1)
      IMD2         = IMD
C--
      BRW2 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      WAT  = WAT/BRW2
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( Z(6).LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*Z(6)/HLM + 1
      ELSE
         ICMB = 2
         JCMB = NHEL*(Z(6)-HLM)/(1-HLM) + 1
      ENDIF
      WAT = WAT*NHEL
      JCMB = MIN(JCMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
C--
C  Select kinematic variables.
C--
      AM(1)   = AMF(IDPTCL(1,1),IDPTCL(2,1),IDPTCL(3,1))
      AM(2)   = AMF(IDPTCL(1,2),IDPTCL(2,2),IDPTCL(3,2))
      AM(3)   = SZM(1)
      AM(4)   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM(5)   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM(6)   = SZM(1)
      AM(7)   = AMF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AM(8)   = AMF(IDPTCL(1,8),IDPTCL(2,8),IDPTCL(3,8))
C--
      AMR(1,1,1) = AMSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4)) 
      AMR(2,1,1) = GTSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AMR(3,1,1) = AMR(1,1,1) + AM(5)
      AMR(1,2,1) = AMSF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AMR(2,2,1) = GTSF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AMR(3,2,1) = AMR(1,2,1) + AM(4)
      AMR(1,3,1) = AMW
      AMR(2,3,1) = GMWTOT
      AMR(3,3,1) = AMR(1,3,1) + AM(3)
      AMR(1,1,2) = AMSF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7)) 
      AMR(2,1,2) = GTSF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AMR(3,1,2) = AMR(1,1,2) + AM(8)
      AMR(1,2,2) = AMSF(IDPTCL(1,8),IDPTCL(2,8),IDPTCL(3,8))
      AMR(2,2,2) = GTSF(IDPTCL(1,8),IDPTCL(2,8),IDPTCL(3,8))
      AMR(3,2,2) = AMR(1,2,2) + AM(7)
      AMR(1,3,2) = AMW
      AMR(2,3,2) = GMWTOT
      AMR(3,3,2) = AMR(1,3,2) + AM(6)
C--
C  Set phase space variables.
C--
      CALL GETXPH(RS,AM,AMX,GMX,AMR,Z(7),IPV,XPHASE,WT)
      IF ( WT.LE.0.D0 ) THEN
         FUNC = 0
         RETURN
      ENDIF
      WAT = WAT*WT
C--
C  Calculate differential cross sections.
C--
      CALL SGCXXF(IDPTCL,IHLCMB,AM,IPV,XPHASE,SG,PV)
      SG     = SG*QED
      FUNC   = SG*WAT
      DSGMDX = SG
      IF ( FUNC.LE.0.D0 ) RETURN
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
      IF ( ABS(QV(3,1)).GT.1.E-3 ) THEN
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
      CALL UVCOPY(4*NP,PVOUT(0,1),PV(0,1))
C--
      EBM   = RS/2
      S     = RS*RS
      CALL PSUMxx(PV(0,3),PV(0,6),QV(0,1))
      PTOP  = SQRT( QV(1,1)**2 + QV(2,1)**2 )
C--
      CALL PSUMxx(PV(0,4),PV(0,5),QV(0,2))
      CALL PSUMxx(QV(0,2),PV(0,3),QV(0,11))
      CALL PSUMxx(PV(0,7),PV(0,8),QV(0,3))
      CALL PSUMxx(QV(0,3),PV(0,6),QV(0,12))
      CALL PMIRRx(QV(0,11),QV(0,11))
      CALL PMIRRx(QV(0,12),QV(0,12))
      CALL BOOSTx(PV(0,5),QV(0,11),QV(0,15))
      CALL BOOSTx(PV(0,8),QV(0,12),QV(0,18))
      CS15 = UDOT3(QV(1,15),QV(1,11))/UABSV3(QV(1,15))/UABSV3(QV(1,11))
      CS18 = UDOT3(QV(1,18),QV(1,12))/UABSV3(QV(1,18))/UABSV3(QV(1,12))
C--
      CALL PSUMxx(PV(0,3),PV(0,5),QV(0,13))
      CALL PSUMxx(PV(0,6),PV(0,8),QV(0,16))
      CALL PMIRRx(QV(0,13),QV(0,13))
      CALL PMIRRx(QV(0,16),QV(0,16))
      CALL BOOSTx(PV(0,5),QV(0,13),QV(0,15))
      CALL BOOSTx(PV(0,8),QV(0,16),QV(0,18))
      CS25 = UDOT3(QV(1,15),QV(1,13))/UABSV3(QV(1,15))/UABSV3(QV(1,13))
      CS28 = UDOT3(QV(1,18),QV(1,16))/UABSV3(QV(1,18))/UABSV3(QV(1,16))
C--      
      COSL  =  PV(3,8)/UABSV3(PV(1,8))
      COSLB = -PV(3,5)/UABSV3(PV(1,5))
      COSJ  =  QV(3,3)/UABSV3(QV(1,3))
      COSJB = -QV(3,2)/UABSV3(QV(1,2))
      PT45  = MAX(SQRT( QV(1,2)**2 + QV(2,2)**2 ),.1)
      PT78  = MAX(SQRT( QV(1,3)**2 + QV(2,3)**2 ),.1)
      ACOP  = ( QV(1,2)*QV(1,3) + QV(2,2)*QV(2,3) )/PT45/PT78
      IF ( ABS(ACOP).GE.1.0 ) ACOP = SIGN(1.0,ACOP)
      ACOP  = xPI - ACOS(ACOP)
      COST  = XPHASE(2,1)
      PHIT  = XPHASE(3,1)
      COSBB = XPHASE(2,2)
      PHIBB = XPHASE(3,2)
      COSB  = XPHASE(2,4)
      PHIB  = XPHASE(3,4)
      COMB  = JCMB + NHEL*(ICMB-1) + .1
      QXP2  = XPHASE(1,2)
      Q42   = XPHASE(1,3)
      Q52   = XPHASE(2,3)
      QXM2  = XPHASE(1,4)
      Q72   = XPHASE(1,5)
      Q82   = XPHASE(2,5)
C-- 
      CALL XHFILL( 'h01', DBLE(COST)        ,FUNC )
      CALL XHFILL( 'h02', DBLE(PHIT)*xRD2DG ,FUNC )
      CALL XHFILL( 'h03', DBLE(COSBB)       ,FUNC )
      CALL XHFILL( 'h04', DBLE(PHIBB)*xRD2DG,FUNC )
      CALL XHFILL( 'h05', DBLE(COSB )       ,FUNC )
      CALL XHFILL( 'h06', DBLE(PHIB )*xRD2DG,FUNC )
      CALL XHFILL( 'h07', DBLE(SQRT(Q42))   ,FUNC )
      CALL XHFILL( 'h08', DBLE(SQRT(Q52))   ,FUNC )
      CALL XHFILL( 'h09', DBLE(PTOP)        ,FUNC )
      CALL XHFILL( 'h10', DBLE(ACOP)*xRD2DG ,FUNC )
      CALL XHFILL( 'h11', DBLE(SQRT(Q72))   ,FUNC )
      CALL XHFILL( 'h12', DBLE(SQRT(Q82))   ,FUNC )
      CALL XHFILL( 'h13', DBLE(COMB)        ,FUNC )
      IF ( IDPTCL(3,5).EQ.1 .AND. IDPTCL(1,5).LE.2 ) THEN
         CALL XHFILL( 'h14', DBLE(PV(0,5))     ,FUNC )
         CALL XHFILL( 'h15', DBLE(COSLB)       ,FUNC )
         CALL XHFILL( 'h26', DBLE(CS25)        ,FUNC )
      ENDIF
      IF ( IDPTCL(3,8).EQ.1 .AND. IDPTCL(1,8).LE.2 ) THEN
         CALL XHFILL( 'h14', DBLE(PV(0,8))     ,FUNC )
         CALL XHFILL( 'h15', DBLE(COSL)        ,FUNC )
         CALL XHFILL( 'h27', DBLE(CS28)        ,FUNC )
      ENDIF
      IF ( IDPTCL(3,5).EQ.2 ) THEN
         CALL XHFILL( 'h16', DBLE(QV(0,2))     ,FUNC )
         CALL XHFILL( 'h17', DBLE(COSJB)       ,FUNC )
      ENDIF
      IF ( IDPTCL(3,8).EQ.2 ) THEN
         CALL XHFILL( 'h16', DBLE(QV(0,3))     ,FUNC )
         CALL XHFILL( 'h17', DBLE(COSJ)        ,FUNC )
      ENDIF
      EVTYP = IDPTCL(3,4) + IDPTCL(3,7) - 0.9
      CALL XHFILL( 'h18', DBLE(EVTYP)          ,FUNC )
      CALL XHFILL( 'h19', DBLE(IMD1+.1)        ,FUNC )
      CALL XHFILL( 'h20', DBLE(IMD2+.1)        ,FUNC )
      CALL XHFILL( 'h21', DBLE(RS/ROOTS)       ,FUNC )
      CALL XHFILL( 'h22', DBLE(CS15)           ,FUNC )
      CALL XHFILL( 'h23', DBLE(CS18)           ,FUNC )
      CALL XHFILL( 'h24', DBLE(SQRT(QXP2))     ,FUNC )
      CALL XHFILL( 'h25', DBLE(SQRT(QXM2))     ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
