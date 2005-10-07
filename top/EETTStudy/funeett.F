C* (Update Record)
C*   95/04/26  K.Fujii	FUNC for e+ e- -> e+ e- t tbar full amplitudes.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
#include "EETTCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 10, NHEL = 8 )
      REAL   *8   SG, PV(0:3,NP)
      REAL   *8   EBM, EMBM, EPBM, EMBM0, EPBM0, PEM, PEP
C--
      REAL   *8   QV(0:3,20)
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,-1,   -1, -1,   -1,+1,-1, +1,-1,+1,
     .             -1,-1,   -1, +1,   -1,+1,-1, +1,-1,+1,
     .             -1,-1,   +1, -1,   -1,+1,-1, +1,-1,+1,
     .             -1,-1,   +1, +1,   -1,+1,-1, +1,-1,+1,
     .             -1,+1,   -1, -1,   -1,+1,-1, +1,-1,+1,
     .             -1,+1,   -1, +1,   -1,+1,-1, +1,-1,+1,
     .             -1,+1,   +1, -1,   -1,+1,-1, +1,-1,+1,
     .             -1,+1,   +1, +1,   -1,+1,-1, +1,-1,+1/
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, -1,   -1,+1,-1, +1,-1,+1,
     .             +1,-1,   -1, +1,   -1,+1,-1, +1,-1,+1,
     .             +1,-1,   +1, -1,   -1,+1,-1, +1,-1,+1,
     .             +1,-1,   +1, +1,   -1,+1,-1, +1,-1,+1,
     .             +1,+1,   -1, -1,   -1,+1,-1, +1,-1,+1,
     .             +1,+1,   -1, +1,   -1,+1,-1, +1,-1,+1,
     .             +1,+1,   +1, -1,   -1,+1,-1, +1,-1,+1,
     .             +1,+1,   +1, +1,   -1,+1,-1, +1,-1,+1/
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
C--
         AMT   = AMF(3,1,2)
         GMT   = GMF(3,1,2)
         AMGT  = AMT*GMT
         AMT2  = AMT*AMT
         AMGT2 = AMGT*AMGT
C--
         GMW   = GMWTOT
         AMGW  = AMW*GMW
         AMW2  = AMW*AMW
         AMGW2 = AMGW*AMGW
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895E0
C--
         XW1MN = BRVDK(IMD1MN-1,1)
         XW1MX = BRVDK(IMD1MX  ,1)
         XW2MN = BRVDK(IMD2MN-1,1)
         XW2MX = BRVDK(IMD2MX  ,1)
      ENDIF
C--
C  Shuffle integration variables.
C--
      DO 10 I = 1, NZZ
         Z(ISHUFL(I)) = ZZ(I)
10    CONTINUE 
C--
C  Set independent variables.
C     Z( 1) : e- beam
C      ( 2) : e+ beam
C      ( 3) : bremsstrahlng
C      ( 4) : helicity
C      ( 5) : xi   
C      ( 6) : eta_1
C      ( 7) : eta_2
C      ( 8) : m(ttbar)**2
C      ( 9) : m(t)**2
C      (10) : m(tbar)**2
C      (11) : m(W+)**2
C      (12) : m(W-)**2
C      (13) : phi_1
C      (14) : phi_2 - phi_1
C      (15) : cos_t        in tt frame
C      (16) : phi_t        in tt frame
C      (17) : cos_b        in t  frame
C      (18) : phi_b   	   in t  frame
C      (19) : cos_bb       in tb frame
C      (20) : phi_bb       in tb frame
C      (21) : cos(theta_f_bar) in W+ rest frame
C      (22) : phi_f_bar        in W+ rest frame
C      (23) : cos(theta_f)     in W- rest frame
C      (24) : phi_f            in W- rest frame
C      (25) : final state combination.
C      (26) : E_beam spread of E-
C      (27) : E_beam spread of E+
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
         CALL PROBEB(Z(26),Z(1),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(Z(27),Z(2),EBM,EPBM,EPBM0,PEP)
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
         IZE = Z(13)*1.D8
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
C     t --> b + W 
C--
      IDPTCL(1, 3) = 1
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 1
      IDPTCL(2, 4) = 2
      IDPTCL(3, 4) = 1
C--
      IDPTCL(1, 5) = 3
      IDPTCL(2, 5) = 2
      IDPTCL(3, 5) = 2
      IDPTCL(1, 8) = 3
      IDPTCL(2, 8) = 2
      IDPTCL(3, 8) = 2
C
C-- W-.
      XW1 = XW1MN + Z(25)*(XW1MX-XW1MN)
      DO 100 IMD = IMD1MN, IMD1MX
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 110
100   CONTINUE
C--
110   IDPTCL(1,6) = MDVDK(2,IMD,1)
      IDPTCL(2,6) = 2
      IDPTCL(3,6) = MDVDK(3,IMD,1)
      IDPTCL(1,7) = MDVDK(1,IMD,1)
      IDPTCL(2,7) = 1
      IDPTCL(3,7) = MDVDK(3,IMD,1)
C--
      BRW1 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD1 = IMD
      WAT  = WAT*(XW1MX-XW1MN)/BRW1
C-- W+.
      XW2 = XW2MN + ((XW1-BRVDK(IMD-1,1))/BRW1)*(XW2MX-XW2MN)
      DO 200 IMD = IMD2MN, IMD2MX
         IF ( XW2.LT.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1, 9) = MDVDK(2,IMD,1)
      IDPTCL(2, 9) = 2
      IDPTCL(3, 9) = MDVDK(3,IMD,1)
      IDPTCL(1,10) = MDVDK(1,IMD,1)
      IDPTCL(2,10) = 1
      IDPTCL(3,10) = MDVDK(3,IMD,1)
C--
      BRW2 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD2 = IMD
      WAT  = WAT*(XW2MX-XW2MN)/BRW2
C--
C  Select helicity combination.
C--
      XHL  = Z(4)
      HLM  = (1-POLE)/2
      IF ( XHL.LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*XHL/HLM + 1
      ELSE
         ICMB = 2
         JCMB = NHEL*(XHL-HLM)/(1-HLM) + 1
      ENDIF
      WAT  = WAT*NHEL
      JCMB = MIN(JCMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
      XHCMB = NHEL*(ICMB-1) + JCMB
C--
C  Return if RS is too small.
C--
      AM3   = AMF(IDPTCL(1, 3),IDPTCL(2, 3),IDPTCL(3, 3))
      AM4   = AMF(IDPTCL(1, 4),IDPTCL(2, 4),IDPTCL(3, 4))
      AM5   = AMF(IDPTCL(1, 5),IDPTCL(2, 5),IDPTCL(3, 5))
      AM6   = AMF(IDPTCL(1, 6),IDPTCL(2, 6),IDPTCL(3, 6))
      AM7   = AMF(IDPTCL(1, 7),IDPTCL(2, 7),IDPTCL(3, 7))
      AM8   = AMF(IDPTCL(1, 8),IDPTCL(2, 8),IDPTCL(3, 8))
      AM9   = AMF(IDPTCL(1, 9),IDPTCL(2, 9),IDPTCL(3, 9))
      AM10  = AMF(IDPTCL(1,10),IDPTCL(2,10),IDPTCL(3,10))
C--
      IF ( RS.LE.AM3+AM4+AM5+AM6+AM7+AM8+AM9+AM10 ) THEN
         FUNC = 0
         RETURN
      ENDIF
C--
C  Set kinematics.
C--
C-- q^2_tt.
      Q32MN = ( AM5 + AM6 + AM7 + AM8 + AM9 + AM10 )**2
      Q32MX = ( RS - ( AM3 + AM4 ) )**2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(8)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
C-- q^2_t.
      Q4MN  = AM5 + AM6 + AM7
      Q4MX  = Q3  - ( AM8 + AM9 + AM10 )
      TH4MN = ATAN((Q4MN-AMT)*(Q4MN+AMT)/AMGT)
      TH4MX = ATAN((Q4MX-AMT)*(Q4MX+AMT)/AMGT)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(9)
      Q42   = AMGT*TAN(TH4) + AMT2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMT2)**2 + AMGT2 )/AMGT
C-- q^2_tb.
      Q5MN  = AM8 + AM9 + AM10
      Q5MX  = Q3  - Q4
      TH5MN = ATAN((Q5MN-AMT)*(Q5MN+AMT)/AMGT)
      TH5MX = ATAN((Q5MX-AMT)*(Q5MX+AMT)/AMGT)
      TH5   = TH5MN + (TH5MX-TH5MN)*Z(10)
      Q52   = AMGT*TAN(TH5) + AMT2
      Q5    = SQRT(Q52)
      WAT   = WAT*(TH5MX-TH5MN)*( (Q52-AMT2)**2 + AMGT2 )/AMGT
C-- q^2_W+.
      Q6MN  = AM6 + AM7
      Q6MX  = Q4  - AM5
      TH6MN = ATAN((Q6MN-AMW)*(Q6MN+AMW)/AMGW)
      TH6MX = ATAN((Q6MX-AMW)*(Q6MX+AMW)/AMGW)
      TH6   = TH6MN + (TH6MX-TH6MN)*Z(11)
      Q62   = AMGW*TAN(TH6) + AMW2
      Q6    = SQRT(Q62)
      WAT   = WAT*(TH6MX-TH6MN)*( (Q62-AMW2)**2 + AMGW2 )/AMGW
C-- q^2_W-.
      Q7MN  = AM9 + AM10
      Q7MX  = Q5  - AM8
      TH7MN = ATAN((Q7MN-AMW)*(Q7MN+AMW)/AMGW)
      TH7MX = ATAN((Q7MX-AMW)*(Q7MX+AMW)/AMGW)
      TH7   = TH7MN + (TH7MX-TH7MN)*Z(12)
      Q72   = AMGW*TAN(TH7) + AMW2
      Q7    = SQRT(Q72)
      WAT   = WAT*(TH7MX-TH7MN)*( (Q72-AMW2)**2 + AMGW2 )/AMGW
C-- Xi_lo.
      XILO = LOG(Q3*(Q3+AM3+AM4)/S)
      XIHI = LOG(1-2*AM3/RS)
      XI    = XILO + (XIHI-XILO)*Z(5)
      RXI   = EXP(XI)
      WAT   = WAT*(XIHI-XILO)
C--
      DM1   = (AME*AME/S)*RXI*RXI/(1-RXI)
      DP1   = (AME*AME/S)
      DM2   = DP1
      DP2   = DM1
C-- Eta_1.
      ETLO  = -LOG( (1+DM1)/DP1 )/2
      ETHI  =  LOG( (1+DP1)/DM1 )/2
      ET1   = ETLO + (ETHI-ETLO)*Z(6)
      WAT   = WAT*(ETHI-ETLO)
C-- Eta_2.
      ETLO  = -LOG( (1+DM2)/DP2 )/2
      ETHI  =  LOG( (1+DP2)/DM2 )/2
      ET2   = ETLO + (ETHI-ETLO)*Z(7)
      WAT   = WAT*(ETHI-ETLO)
C-- Phi_1.
      FI1   = Z(13)*x2PI
C-- Phi_2.
      FI2   = MOD(Z(14)*x2PI + FI1,x2PI)
C--
C  Set XPHASE.
C--
      XPHASE(1,1) = S
      XPHASE(2,1) = ET1
      XPHASE(3,1) = FI1
      XPHASE(1,2) = XI
      XPHASE(2,2) = ET2
      XPHASE(3,2) = FI2
      XPHASE(1,3) = Q32
      XPHASE(2,3) = -1 + 2*Z(15)
      XPHASE(3,3) = Z(16)*x2PI
      XPHASE(1,4) = Q42
      XPHASE(2,4) = -1 + 2*Z(17)
      XPHASE(3,4) = Z(18)*x2PI
      XPHASE(1,5) = Q52
      XPHASE(2,5) = -1 + 2*Z(19)
      XPHASE(3,5) = Z(20)*x2PI
      XPHASE(1,6) = Q62
      XPHASE(2,6) = -1 + 2*Z(21)
      XPHASE(3,6) = Z(22)*x2PI
      XPHASE(1,7) = Q72
      XPHASE(2,7) = -1 + 2*Z(23)
      XPHASE(3,7) = Z(24)*x2PI
C--
      WAT   = WAT*2**(NP-5)*(x2PI)**(NP-3)
C--
C  Calculate differential cross sections.
C--
      CALL SGEETTF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
         CALL UVCOPY(8*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
      RDDG   = 180/xPI
C--
      EXP2ET  = EXP(2*ET1)
      SH1     = MIN( SQRT( (1+DP1+DM1)/(1+  EXP2ET) - DM1 ), 1.D0 )
      CH1     = MIN( SQRT( (1+DP1+DM1)/(1+1/EXP2ET) - DP1 ), 1.D0 )
      SN1     = 2*SH1*CH1
      CS1     = (1+DP1+DM1)*TANH(ET1) - DP1 + DM1
C--
      EXP2ET  = EXP(2*ET2)
      SH2     = MIN( SQRT( (1+DP2+DM2)/(1+  EXP2ET) - DM2 ), 1.D0 )
      CH2     = MIN( SQRT( (1+DP2+DM2)/(1+1/EXP2ET) - DP2 ), 1.D0 )
      SN2     = 2*SH2*CH2
      CS2     = (1+DP2+DM2)*TANH(ET2) - DP2 + DM2
C--
      CS12    = CS1*CS2 + SN1*SN2*COS(FI2-FI1)
      IF ( ET1.GT.-ET2 ) THEN
         X1      = 1 - RXI
         X2      = ( RXI - Q32/S )/( 1 - X1*(1-CS12)/2 )
      ELSE
         X2      = 1 - RXI
         X1      = ( RXI - Q32/S )/( 1 - X2*(1-CS12)/2 )
      ENDIF
C--
      CALL PSUBxx(PV(0,3),PV(0,1),QV(0,2))
      AQV     = SQRT(QV(1,2)**2+QV(2,2)**2+QV(3,2)**2)
      TEE     = (QV(0,2)-AQV)*(QV(0,2)+AQV)
C--
      CALL PSUMxx(PV(0,6),PV(0,7),QV(0,4))
      CALL PSUMxx(QV(0,4),PV(0,5),QV(0,4))
      XTT1    = 2*SQRT(QV(1,4)**2+QV(2,4)**2)/ROOTS
C--
      CALL PSUMxx(PV(0,9),PV(0,10),QV(0,5))
      CALL PSUMxx(QV(0,5),PV(0, 8),QV(0,5))
      XTT2    = 2*SQRT(QV(1,5)**2+QV(2,5)**2)/ROOTS
C--
      XNU     = 2*PV(0,3)/ROOTS
      XNB     = 2*PV(0,4)/ROOTS
      XTP     = 2*QV(0,4)/ROOTS
      XTB     = 2*QV(0,5)/ROOTS
      CSTP    = XPHASE(2,3)
      FITP    = XPHASE(3,3)
      CSFD    = XPHASE(2,4)
      FIFD    = XPHASE(3,4)
      CSFDB   = XPHASE(2,5)
      FIFDB   = XPHASE(3,5)
      DFI     = MOD(FI1-FI2+x2PI,x2PI)
C--
      CALL XHFILL( 'h01', DBLE(X1)          ,FUNC )
      CALL XHFILL( 'h02', DBLE(CS1)         ,FUNC )
      CALL XHFILL( 'h03', DBLE(FI1)*RDDG    ,FUNC )
      CALL XHFILL( 'h04', DBLE(X2)          ,FUNC )
      CALL XHFILL( 'h05', DBLE(CS2)         ,FUNC )
      CALL XHFILL( 'h06', DBLE(DFI)*RDDG    ,FUNC )
      CALL XHFILL( 'h07', DBLE(Q3/ROOTS)    ,FUNC )
      CALL XHFILL( 'h08', DBLE(CSTP)        ,FUNC )
      CALL XHFILL( 'h09', DBLE(FITP)*RDDG   ,FUNC )
      CALL XHFILL( 'h10', DBLE(Q4)          ,FUNC )
      CALL XHFILL( 'h11', DBLE(CSFD)        ,FUNC )
      CALL XHFILL( 'h12', DBLE(FIFD)*RDDG   ,FUNC )
      CALL XHFILL( 'h13', DBLE(Q5)          ,FUNC )
      CALL XHFILL( 'h14', DBLE(CSFDB)       ,FUNC )
      CALL XHFILL( 'h15', DBLE(FIFDB)*RDDG  ,FUNC )
      CALL XHFILL( 'h16', DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 'h17', DBLE(IMD1+  .1)   ,FUNC )
      CALL XHFILL( 'h17', DBLE(IMD2+12.1)   ,FUNC )
      CALL XHFILL( 'h18', DBLE(XHCMB+.1)    ,FUNC )
      CALL XHFILL( 'h19', DBLE(ET1)         ,FUNC )
      CALL XHFILL( 'h20', DBLE(ET2)         ,FUNC )
      CALL XHFILL( 'h21', DBLE(XTT1)        ,FUNC )
      CALL XHFILL( 'h22', DBLE(XTT2)        ,FUNC )
      CALL XHFILL( 'h23', DBLE(XNU)         ,FUNC )
      CALL XHFILL( 'h24', DBLE(XNB)         ,FUNC )
      CALL XHFILL( 'h25', DBLE(XTP)         ,FUNC )
      CALL XHFILL( 'h26', DBLE(XTB)         ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
