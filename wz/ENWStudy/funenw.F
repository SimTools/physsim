C* (Update Record)
C*   95/04/22  K.Fujii	FUNC for e+ e- --> e nu W full amplitudes.
C*
      REAL FUNCTION FUNC*8(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "ENWCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6, NHEL = 4 )
      REAL   *8   SG, PV(0:3,NP)
C--
      REAL   *8   QV(0:3,20)
      REAL   *8   Q3, Q32, Q3MN, Q3MX, TH3MN, TH3MX, TH3
      REAL   *8   WAT
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,-1,   -1, +1,   +1,-1,
     .             -1,-1,   +1, +1,   +1,-1,
     .             -1,+1,   -1, +1,   +1,-1,
     .             -1,+1,   +1, +1,   +1,-1/
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, +1,   +1,-1,
     .             +1,-1,   +1, +1,   +1,-1,
     .             +1,+1,   -1, +1,   +1,-1,
     .             +1,+1,   +1, +1,   +1,-1/
C--
      INTEGER*4  JHEL(NP,NHEL,2)
      DATA ( ( JHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             +1,+1,   -1, +1,   +1,-1,
     .             +1,+1,   +1, +1,   +1,-1,
     .             -1,+1,   -1, +1,   +1,-1,
     .             -1,+1,   +1, +1,   +1,-1/
      DATA ( ( JHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, +1,   +1,-1,
     .             +1,-1,   +1, +1,   +1,-1,
     .             -1,-1,   -1, +1,   +1,-1,
     .             -1,-1,   +1, +1,   +1,-1/
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
C--
         GMW   = GMWTOT
         AMGW  = AMW*GMW
         AMW2  = AMW*AMW
         AMGW2 = AMGW*AMGW
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895D0
C--
         XW1MN = BRVDK(IMDMN-1,1)
         XW1MX = BRVDK(IMDMX  ,1)
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
C      ( 4) : helicity
C      ( 5) : xi   (E_e)
C      ( 6) : m(W)**2
C      ( 7) : eta  (cos_e)
C      ( 8) : zeta (cos_n) in (nu W) frame
C      ( 9) : phi_e
C      (10) : phi_n   	   in (nu W) frame
C      (11) : cos_fd       in W+     frame
C      (12) : phi_fd       in W+     frame
C      (13) : final state combination.
C      (14) : E_beam spread of E-
C      (15) : E_beam spread of E+
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Beam energy.
C--
      RS        = ROOTS
      EMBM      = RS/2
      EPBM      = RS/2
      QED       = 1
C--
C  First decide e- and e+ energy.
C--
      IF ( ISRBM.EQ.3 ) THEN
         EBM = ROOTS/2
         CALL PROBEB(Z(14),Z(1),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(Z(15),Z(2),EBM,EPBM,EPBM0,PEP)
         RS  = 2*SQRT(EMBM*EPBM)
         IF ( RS.LT.5.D0 ) THEN
            FUNC = 0
            RETURN
         ENDIF
      ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(14)*1.D4
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
C--      
      S      = RS*RS
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
      IDPTCL(1, 3) = 1
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 1
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = 1
C-- W.
      XW1 = XW1MN + Z(13)*(XW1MX-XW1MN)
      DO 200 IMD = IMDMN, IMDMX
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,5) = MDVDK(2,IMD,1)
      IDPTCL(2,5) = 2
      IDPTCL(3,5) = MDVDK(3,IMD,1)
      IDPTCL(1,6) = MDVDK(1,IMD,1)
      IDPTCL(2,6) = 1
      IDPTCL(3,6) = MDVDK(3,IMD,1)
C--
      BRW1 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      WAT  = WAT*(XW1MX-XW1MN)/BRW1
C--
C  Select CP.
C--
      IF ( Z(4).LT.0.5D0 ) THEN
         NCPFLG = +1
         XHL    = 2*Z(4)
      ELSE
         NCPFLG = -1
         XHL    = 2*Z(4) - 1
      ENDIF
      WAT   = WAT*2
C>>>
C      XHL    = Z(4)
C      WAT    = WAT/2
C      NCPFLG = +1
C>>>
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( XHL.LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*XHL/HLM + 1
      ELSE
         ICMB = 2
         JCMB = NHEL*(XHL-HLM)/(1-HLM) + 1
      ENDIF
      WAT  = WAT*NHEL
      JCMB = MIN(JCMB,NHEL)
      IF ( NCPFLG.GT.0 ) THEN
         CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
      ELSE
         CALL UVCOPY(NP,JHEL(1,JCMB,ICMB),IHLCMB(1))
      ENDIF
      XHCMB = NHEL*(ICMB-1) + JCMB + NHEL*(1-NCPFLG)
C--
C  Set kinematics.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
C-- q^2_W.
      Q3MN  = AM5 + AM6
      Q3MX  = RS - ( AM3 + AM4 )
      TH3MN = ATAN((Q3MN-AMW)*(Q3MN+AMW)/AMGW)
      TH3MX = ATAN((Q3MX-AMW)*(Q3MX+AMW)/AMGW)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(6)
      Q32   = AMGW*TAN(TH3) + AMW2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMW2)**2 + AMGW2 )/AMGW
C-- Xi_lo.
      XILO  = LOG( (Q3-AME)*(Q3+AME)/S )
      XIHI  = LOG( 1 - 2*AME/RS )
      XI    = XILO + (XIHI-XILO)*Z(5)
      WAT   = WAT*(XIHI-XILO)
      EXPXI = EXP(XI)
C-- Eta.
      DME   = (AME*AME/S)*EXPXI**2/(1-EXPXI)
      DPE   = 1
      ETLO  = -LOG( (1+DME)/DPE )/2
      ETHI  =  LOG( (1+DPE)/DME )/2
      ETA   = ETLO + (ETHI-ETLO)*Z(7)
      WAT   = WAT*(ETHI-ETLO)
C--
C  Set XPHASE.
C--
      XPHASE(1,1) = S
      XPHASE(2,1) = ETA
      XPHASE(3,1) = Z(9)*x2PI
      XPHASE(1,2) = XI
      XPHASE(2,2) = Z(8)
      XPHASE(3,2) = Z(10)*x2PI
      XPHASE(1,3) = Q32
      XPHASE(2,3) = -1 + 2*Z(11)
      XPHASE(3,3) = Z(12)*x2PI
C--
      WAT   = WAT*2*(x2PI)**3
C--
C  Calculate differential cross sections.
C--
      CALL SGENWF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG     = SG*QED
      FUNC   = SG*WAT
      IF ( FUNC.LE.0.D0 ) RETURN
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
      IF ( ABS(QV(3,1)).GT.1.D-3 ) THEN
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(2*4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Flip everything when NCPFLG = -1.
C--
      IF ( NCPFLG.EQ.-1 ) THEN
         DO 600 IP = 3, NP
            CALL PMIRRx(PVOUT(0,IP),PVOUT(0,IP))
            IHLCMB(IP) = -IHLCMB(IP)
600      CONTINUE            
      ENDIF
C--
C  Fill hists. and plots.
C--
      RDDG   = 180/xPI
      QWN2   = S*EXPXI + AME*AME
      QWN    = SQRT(QWN2)
      EXP2ET = EXP(2*XPHASE(2,1))
      SH1    = MIN( SQRT( (1+DPE+DME)/(1+EXP2ET) - DME ), 1.D0 )
      DLT    = 2*AMW*AMW/((QWN-AMW)*(QWN+AMW))
      OPDZT  = ( 1 + 1/DLT )**XPHASE(2,2)
      SH2    = MIN( SQRT( DLT*( OPDZT - 1 ) ), 1.D0 )
C--
      XE    = 1 + AME*AME/S - QWN2/S
      COSE  = ( 1 - xSQRT2*SH1 )*( 1 + xSQRT2*SH1 )
      PHIE  = XPHASE(3,1)
      COSN  = ( 1 - xSQRT2*SH2 )*( 1 + xSQRT2*SH2 )
      PHIN  = XPHASE(3,2)
      QW    = Q3
      CSDA  = XPHASE(2,3)
      PHID  = XPHASE(3,3)
C--
      QV(0,2) = PV(0,3) - PV(0,1)
      QV(1,2) = PV(1,3) - PV(1,1)
      QV(2,2) = PV(2,3) - PV(2,1)
      QV(3,2) = PV(3,3) - PV(3,1)
      AQV     = SQRT(QV(1,2)**2+QV(2,2)**2+QV(3,2)**2)
      TEE     = (QV(0,2)-AQV)*(QV(0,2)+AQV)
C--
      QV(0,5) = PV(0,5) + PV(0,6)
      QV(1,5) = PV(1,5) + PV(1,6)
      QV(2,5) = PV(2,5) + PV(2,6)
      QV(3,5) = PV(3,5) + PV(3,6)
      XTW     = 2*SQRT(QV(1,5)**2+QV(2,5)**2)/ROOTS
C--
      XEL     = 2*PV(0,3)/ROOTS
      XNU     = 2*PV(0,4)/ROOTS
      XWP     = 2*QV(0,5)/ROOTS
C--
      CALL XHFILL(  1, COSE              ,FUNC )
      CALL XHFILL(  2, PHIE*RDDG         ,FUNC )
      CALL XHFILL(  3, XE                ,FUNC )
      CALL XHFILL(  4, QW                ,FUNC )
      CALL XHFILL(  5, COSN              ,FUNC )
      CALL XHFILL(  6, PHIN*RDDG         ,FUNC )
      CALL XHFILL(  7, CSDA              ,FUNC )
      CALL XHFILL(  8, PHID*RDDG         ,FUNC )
      CALL XHFILL(  9, RS/ROOTS          ,FUNC )
      CALL XHFILL( 10, IMD+.1D0          ,FUNC )
      CALL XHFILL( 11, XHCMB+.1D0        ,FUNC )
      CALL XHFILL( 12, ETA               ,FUNC )
      CALL XHFILL( 13, XTW               ,FUNC )
      CALL XHFILL( 14, XEL               ,FUNC )
      CALL XHFILL( 15, XNU               ,FUNC )
      CALL XHFILL( 16, XWP               ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
