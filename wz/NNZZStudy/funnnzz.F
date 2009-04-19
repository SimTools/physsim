C* (Update Record)
C*   95/04/26  K.Fujii	FUNC for e- e+ --> nu nub Z Z full amplitudes.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
#include "NNZZCNF.inc"
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
     .             -1,+1,   -1, +1,   -1,+1, -1,+1,
     .             -1,+1,   -1, +1,   -1,+1, +1,-1,
     .             -1,+1,   -1, +1,   +1,-1, -1,+1,
     .             -1,+1,   -1, +1,   +1,-1, +1,-1/
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, +1,   -1,+1, -1,+1,
     .             +1,-1,   -1, +1,   -1,+1, +1,-1,
     .             +1,-1,   -1, +1,   +1,-1, -1,+1,
     .             +1,-1,   -1, +1,   +1,-1, +1,-1/
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
         GMV   = GMZTOT
         AMV   = AMZ
         AMGV  = AMV*GMV
         AMV2  = AMV*AMV
         AMGV2 = AMGV*AMGV
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895E0
C--
         XZ1MN = BRVDK(IMD1MN-1,2)
         XZ1MX = BRVDK(IMD1MX  ,2)
         XZ2MN = BRVDK(IMD2MN-1,2)
         XZ2MX = BRVDK(IMD2MX  ,2)
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
C      ( 8) : m(ZZ)**2
C      ( 9) : m(Z1)**2
C      (10) : m(Z2)**2
C      (11) : phi_1
C      (12) : phi_2 - phi_1
C      (13) : cos_Z1       in ZZ frame
C      (14) : phi_Z1       in ZZ frame
C      (15) : cos_f        in Z1 frame
C      (16) : phi_f   	   in Z1 frame
C      (17) : cos_f        in Z2 frame
C      (18) : phi_f        in Z2 frame
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
      IDPTCL(2, 3) = 1
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 1
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = 1
C-- Z1.
      XZ1 = XZ1MN + Z(19)*(XZ1MX-XZ1MN)
      DO 100 IMD = IMD1MN, IMD1MX
         IF ( XZ1.LT.BRVDK(IMD,2) )              GO TO 110
100   CONTINUE
C--
110   IDPTCL(1,5) = MDVDK(1,IMD,2)
      IDPTCL(2,5) = MDVDK(2,IMD,2)
      IDPTCL(3,5) = MDVDK(3,IMD,2)
      IDPTCL(1,6) = MDVDK(1,IMD,2)
      IDPTCL(2,6) = MDVDK(2,IMD,2)
      IDPTCL(3,6) = MDVDK(3,IMD,2)
C--
      BRZ1 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD1 = IMD
      WAT  = WAT*(XZ1MX-XZ1MN)/BRZ1
C-- Z2.
      XZ2 = XZ2MN + ((XZ1-BRVDK(IMD-1,2))/BRZ1)*(XZ2MX-XZ2MN)
      DO 200 IMD = IMD2MN, IMD2MX
         IF ( XZ2.LT.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,7) = MDVDK(1,IMD,2)
      IDPTCL(2,7) = MDVDK(2,IMD,2)
      IDPTCL(3,7) = MDVDK(3,IMD,2)
      IDPTCL(1,8) = MDVDK(1,IMD,2)
      IDPTCL(2,8) = MDVDK(2,IMD,2)
      IDPTCL(3,8) = MDVDK(3,IMD,2)
C--
      BRZ2 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD2 = IMD
      WAT  = WAT*(XZ2MX-XZ2MN)/BRZ2
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
C-- q^2_ZZ
      Q32MN = ( AM5 + AM6 + AM7 + AM8 )**2
      Q32MX = ( RS - ( AM3 + AM4 ) )**2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(8)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
C-- q^2_Z1.
      Q4MN  = AM5 + AM6
      Q4MX  = Q3  - ( AM7 + AM8 )
      TH4MN = ATAN((Q4MN-AMV)*(Q4MN+AMV)/AMGV)
      TH4MX = ATAN((Q4MX-AMV)*(Q4MX+AMV)/AMGV)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(9)
      Q42   = AMGV*TAN(TH4) + AMV2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMV2)**2 + AMGV2 )/AMGV
C-- q^2_Z2.
      Q5MN  = AM7 + AM8
      Q5MX  = Q3  - Q4
      TH5MN = ATAN((Q5MN-AMV)*(Q5MN+AMV)/AMGV)
      TH5MX = ATAN((Q5MX-AMV)*(Q5MX+AMV)/AMGV)
      TH5   = TH5MN + (TH5MX-TH5MN)*Z(10)
      Q52   = AMGV*TAN(TH5) + AMV2
      Q5    = SQRT(Q52)
      WAT   = WAT*(TH5MX-TH5MN)*( (Q52-AMV2)**2 + AMGV2 )/AMGV
C-- Xi_lo.
      XILO = LOG(Q3*(Q3+AM3+AM4)/S)
      XIHI = LOG(1-2*AM3/RS)
      XI    = XILO + (XIHI-XILO)*Z(5)
      RXI   = EXP(XI)
      WAT   = WAT*(XIHI-XILO)
C--
C For gamma gamma
C      DM1   = (AME*AME/S)*RXI*RXI/(1-RXI)
C      DP1   = (AME*AME/S)
C For W W fusion
C
      DM1   = AMW*AMW/S
      DP1   = 1
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
      CALL SGNNZZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
      CALL PSUMxx(PV(0,5),PV(0,6),QV(0,4))
      XTW1    = 2*SQRT(QV(1,4)**2+QV(2,4)**2)/ROOTS
C--
      CALL PSUMxx(PV(0,7),PV(0,8),QV(0,5))
      XTW2    = 2*SQRT(QV(1,5)**2+QV(2,5)**2)/ROOTS
C--
      XEM     = 2*PV(0,3)/ROOTS
      XEP     = 2*PV(0,4)/ROOTS
      XWM     = 2*QV(0,4)/ROOTS
      XWP     = 2*QV(0,5)/ROOTS
      CSWM    = XPHASE(2,3)
      FIWM    = XPHASE(3,3)
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
      CALL XHFILL( 'h08', DBLE(CSWM)        ,FUNC )
      CALL XHFILL( 'h09', DBLE(FIWM)*RDDG   ,FUNC )
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
      CALL XHFILL( 'h21', DBLE(XTW1)        ,FUNC )
      CALL XHFILL( 'h22', DBLE(XTW2)        ,FUNC )
      CALL XHFILL( 'h23', DBLE(XEM)         ,FUNC )
      CALL XHFILL( 'h24', DBLE(XEP)         ,FUNC )
      CALL XHFILL( 'h25', DBLE(XWM)         ,FUNC )
      CALL XHFILL( 'h26', DBLE(XWP)         ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
