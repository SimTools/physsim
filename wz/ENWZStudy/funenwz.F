C* (Update Record)
C*                                              __
C*   2008/11/28  K.Fujii  FUNC for e+ e- --> e- nu_e W+ Z full amplitudes.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
#include "ENWZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 8, NHEL = 4 )
      REAL   *8   SG, PV(0:3,NP)
      REAL   *8   EBM, EMBM, EPBM, EMBM0, EPBM0, PEM, PEP
C--
      REAL   *8   QV(0:3,20)
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,+1,   -1, +1,   +1,-1, +1,-1,
     .             -1,+1,   -1, +1,   +1,-1, -1,+1,
     .             -1,+1,   +1, +1,   +1,-1, +1,-1,
     .             -1,+1,   +1, +1,   +1,-1, -1,+1/
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,+1,   -1, +1,   +1,-1, +1,-1,
     .             +1,+1,   -1, +1,   +1,-1, -1,+1,
     .             +1,+1,   +1, +1,   +1,-1, +1,-1,
     .             +1,+1,   +1, +1,   +1,-1, -1,+1/
C--
      INTEGER*4  JHEL(NP,NHEL,2)
      DATA ( ( JHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,+1,   -1, +1,   +1,-1, +1,-1,
     .             -1,+1,   -1, +1,   +1,-1, -1,+1,
     .             +1,+1,   +1, +1,   +1,-1, +1,-1,
     .             +1,+1,   +1, +1,   +1,-1, -1,+1/
      DATA ( ( JHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             -1,-1,   -1, +1,   +1,-1, +1,-1,
     .             -1,-1,   -1, +1,   +1,-1, -1,+1,
     .             +1,-1,   +1, +1,   +1,-1, +1,-1,
     .             +1,-1,   +1, +1,   +1,-1, -1,+1/
C--
      DATA NCALL /  0 /
      EXTERNAL   XHFILL, DHFILL
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
         GMZ   = GMZTOT
         AMGZ  = AMZ*GMZ
         AMZ2  = AMZ*AMZ
         AMGZ2 = AMGZ*AMGZ
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895E0
C--
         XWPMN = BRVDK(IMD1MN-1,1)
         XWPMX = BRVDK(IMD1MX  ,1)
         XZ0MN = BRVDK(IMD2MN-1,2)
         XZ0MX = BRVDK(IMD2MX  ,2)
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
C      ( 8) : m(WZ)**2
C      ( 9) : m(W+)**2
C      (10) : m(Z)**2
C      (11) : phi_W
C      (12) : phi_W - phi_Z
C      (13) : cos_W+       in WZ frame
C      (14) : phi_W+       in WZ frame
C      (15) : cos_fdb      in W+ frame
C      (16) : phi_fdb 	   in W+ frame
C      (17) : cos_f        in Z  frame
C      (18) : phi_f        in Z  frame
C      (19) : final state combination.
C      (20) : E_beam spread of E-
C      (21) : E_beam spread of E+
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
         call jsfbeamgend(Z(20), Z(21), Z(1), Z(2), ebm,
     >                    embm, epbm, embm0, epbm0, dpdebm )
         RS  = 2*SQRT(EMBM*EPBM)
         IF ( RS.LT.5.D0 ) THEN
            FUNC = 0
            RETURN
         ENDIF
         wat = wat*dpdebm
      ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(20)*1.D8
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
C-- W+.
      XWP = XWPMN + Z(19)*(XWPMX-XWPMN)
      DO 100 IMD = IMD1MN, IMD1MX
         IF ( XWP.LE.BRVDK(IMD,1) )              GO TO 110
100   CONTINUE
C--
110   IDPTCL(1,5) = MDVDK(2,IMD,1)
      IDPTCL(2,5) = 2
      IDPTCL(3,5) = MDVDK(3,IMD,1)
      IDPTCL(1,6) = MDVDK(1,IMD,1)
      IDPTCL(2,6) = 1
      IDPTCL(3,6) = MDVDK(3,IMD,1)
C--
      BRWP = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD1 = IMD
      WAT  = WAT*(XWPMX-XWPMN)/BRWP
C-- W+.
      XZ0 = XZ0MN + ((XWP-BRVDK(IMD-1,1))/BRWP)*(XZ0MX-XZ0MN)
      DO 200 IMD = IMD2MN, IMD2MX
         IF ( XZ0.LE.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,7) = MDVDK(1,IMD,2)
      IDPTCL(2,7) = MDVDK(2,IMD,2)
      IDPTCL(3,7) = MDVDK(3,IMD,2)
      IDPTCL(1,8) = MDVDK(1,IMD,2)
      IDPTCL(2,8) = MDVDK(2,IMD,2)
      IDPTCL(3,8) = MDVDK(3,IMD,2)
C--
      BRZ0 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD2 = IMD
      WAT  = WAT*(XZ0MX-XZ0MN)/BRZ0
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
      IF ( NCPFLG.GT.0 ) THEN
        CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
      ELSE
        CALL UVCOPY(NP,JHEL(1,JCMB,ICMB),IHLCMB(1))
      ENDIF
      XHCMB = NHEL*(ICMB-1) + JCMB + NHEL*(1-NCPFLG)
C--
C  Return if RS is too small.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AM7   = AMF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AM8   = AMF(IDPTCL(1,8),IDPTCL(2,8),IDPTCL(3,8))
C--
      IF ( RS.LE.AM3+AM4+AM5+AM6+AM7+AM8 ) THEN
         FUNC = 0
         RETURN
      ENDIF
C--
C  Set kinematics.
C--
C-- q^2_WZ.
      Q32MN = ( AM5 + AM6 + AM7 + AM8 )**2
      Q32MX = ( RS - ( AM3 + AM4 ) )**2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(8)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
C-- q^2_W+.
      Q4MN  = AM5 + AM6
      Q4MX  = Q3  - ( AM7 + AM8 )
#ifndef __ZEROWIDTH__
      TH4MN = ATAN((Q4MN-AMW)*(Q4MN+AMW)/AMGW)
      TH4MX = ATAN((Q4MX-AMW)*(Q4MX+AMW)/AMGW)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(9)
      Q42   = AMGW*TAN(TH4) + AMW2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMW2)**2 + AMGW2 )/AMGW
#else
      Q42   = AMW**2
      Q4    = AMW
      WAT   = WAT*xPI*AMGW
#endif
C-- q^2_Z0.
      Q5MN  = AM7 + AM8
      Q5MX  = Q3  - Q4
#ifndef __ZEROWIDTH__
      TH5MN = ATAN((Q5MN-AMZ)*(Q5MN+AMZ)/AMGZ)
      TH5MX = ATAN((Q5MX-AMZ)*(Q5MX+AMZ)/AMGZ)
      TH5   = TH5MN + (TH5MX-TH5MN)*Z(10)
      Q52   = AMGZ*TAN(TH5) + AMZ2
      Q5    = SQRT(Q52)
      WAT   = WAT*(TH5MX-TH5MN)*( (Q52-AMZ2)**2 + AMGZ2 )/AMGZ
#else
      Q52   = AMZ**2
      Q5    = AMZ
      WAT   = WAT*xPI*AMGZ
#endif
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
      FI1   = Z(11)*x2PI
C-- Phi_2.
      FI2   = MOD(Z(12)*x2PI + FI1,x2PI)
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
      XPHASE(2,3) = -1 + 2*Z(13)
      XPHASE(3,3) = Z(14)*x2PI
      XPHASE(1,4) = Q42
      XPHASE(2,4) = -1 + 2*Z(15)
      XPHASE(3,4) = Z(16)*x2PI
      XPHASE(1,5) = Q52
      XPHASE(2,5) = -1 + 2*Z(17)
      XPHASE(3,5) = Z(18)*x2PI
C--
      WAT   = WAT*2**3*(x2PI)**5
C--
C  Calculate differential cross sections.
C--
      CALL SGENWZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
C  Flip everything when NCPFLG = -1.
C--
      IF ( NCPFLG.EQ.-1 ) THEN
         IHLCMB(1) = -IHLCMB(1)
         IHLCMB(2) = -IHLCMB(2)
         DO 600 IP = 3, NP
            CALL PMIRRx(PVOUT(0,IP),PVOUT(0,IP))
            IHLCMB(IP) = -IHLCMB(IP)
600      CONTINUE            
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
      CALL PSUMxx(PV(0,5),PV(0,6),QV(0,4))
      XTWP    = 2*SQRT(QV(1,4)**2+QV(2,4)**2)/ROOTS
C--
      CALL PSUMxx(PV(0,7),PV(0,8),QV(0,5))
      XTZ0    = 2*SQRT(QV(1,5)**2+QV(2,5)**2)/ROOTS
C--
      XEM     = 2*PV(0,3)/ROOTS
      XEP     = 2*PV(0,4)/ROOTS
      XWP     = 2*QV(0,4)/ROOTS
      XZ0     = 2*QV(0,5)/ROOTS
      CSWP    = XPHASE(2,3)
      FIWP    = XPHASE(3,3)
      CSFDB   = XPHASE(2,4)
      FIFDB   = XPHASE(3,4)
      CSFD    = XPHASE(2,5)
      FIFD    = XPHASE(3,5)
      DFI     = MOD(FI1-FI2+x2PI,x2PI)
C--
      CALL XHFILL( 'h01', DBLE(X1)          ,FUNC )
      CALL XHFILL( 'h02', DBLE(CS1)         ,FUNC )
      CALL XHFILL( 'h03', DBLE(FI1)*RDDG    ,FUNC )
      CALL XHFILL( 'h04', DBLE(X2)          ,FUNC )
      CALL XHFILL( 'h05', DBLE(CS2)         ,FUNC )
      CALL XHFILL( 'h06', DBLE(DFI)*RDDG    ,FUNC )
      CALL XHFILL( 'h07', DBLE(Q3/ROOTS)    ,FUNC )
      CALL XHFILL( 'h08', DBLE(CSWP)        ,FUNC )
      CALL XHFILL( 'h09', DBLE(FIWP)*RDDG   ,FUNC )
      CALL XHFILL( 'h10', DBLE(Q4)          ,FUNC )
      CALL XHFILL( 'h11', DBLE(CSFDB)       ,FUNC )
      CALL XHFILL( 'h12', DBLE(FIFDB)*RDDG  ,FUNC )
      CALL XHFILL( 'h13', DBLE(Q5)          ,FUNC )
      CALL XHFILL( 'h14', DBLE(CSFD)        ,FUNC )
      CALL XHFILL( 'h15', DBLE(FIFD)*RDDG   ,FUNC )
      CALL XHFILL( 'h16', DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 'h17', DBLE(IMD1+  .1)   ,FUNC )
      CALL XHFILL( 'h17', DBLE(IMD2+12.1)   ,FUNC )
      CALL XHFILL( 'h18', DBLE(XHCMB+.1)    ,FUNC )
      CALL XHFILL( 'h19', DBLE(ET1)         ,FUNC )
      CALL XHFILL( 'h20', DBLE(ET2)         ,FUNC )
      CALL XHFILL( 'h21', DBLE(XTWP)        ,FUNC )
      CALL XHFILL( 'h22', DBLE(XTZ0)        ,FUNC )
      CALL XHFILL( 'h23', DBLE(XEM)         ,FUNC )
      CALL XHFILL( 'h24', DBLE(XEP)         ,FUNC )
      CALL XHFILL( 'h25', DBLE(XWP)         ,FUNC )
      CALL XHFILL( 'h26', DBLE(XZ0)         ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
