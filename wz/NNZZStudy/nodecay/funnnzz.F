C* (Update Record)
C*   2005/03/06  K.Fujii	FUNC for e+ e- --> n nb Z Z full amplitudes.
C*
      REAL FUNCTION FUNC*8(ZZ)
 
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
      PARAMETER  ( NP = 6, NHEL = 9 )
      REAL   *8   SG, PV(0:3,NP)
      REAL   *8   EBM, EMBM, EPBM, EMBM0, EPBM0, PEM, PEP
C--
      REAL   *8   QV(0:3,20)
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,+1,   -1, +1,   -1,-1,
     .             -1,+1,   -1, +1,   -1, 0,
     .             -1,+1,   -1, +1,   -1,+1,
     .             -1,+1,   -1, +1,    0,-1,
     .             -1,+1,   -1, +1,    0, 0,
     .             -1,+1,   -1, +1,    0,+1,
     .             -1,+1,   -1, +1,   +1,-1,
     .             -1,+1,   -1, +1,   +1, 0,
     .             -1,+1,   -1, +1,   +1,+1 /
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, +1,   -1,-1,
     .             +1,-1,   -1, +1,   -1, 0,
     .             +1,-1,   -1, +1,   -1,+1,
     .             +1,-1,   -1, +1,    0,-1,
     .             +1,-1,   -1, +1,    0, 0,
     .             +1,-1,   -1, +1,    0,+1,
     .             +1,-1,   -1, +1,   +1,-1,
     .             +1,-1,   -1, +1,   +1, 0,
     .             +1,-1,   -1, +1,   +1,+1 /
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
         AME   = AMF(1,2,1)
         xALF0 = 1/137.0359895E0
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
C      ( 9) : phi_1
C      (10) : phi_2 - phi_1
C      (11) : cos_Z1        in ZZ frame
C      (12) : phi_Z1        in ZZ frame
C      (13) : E_beam spread of E-
C      (14) : E_beam spread of E+
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
         CALL PROBEB(Z(13),Z(1),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(Z(14),Z(2),EBM,EPBM,EPBM0,PEP)
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
         IZE = Z(9)*1.D8
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
      AM5   = AMZ
      AM6   = AMZ
C--
      IF ( RS.LE.AM3+AM4+AM5+AM6 ) THEN
         FUNC = 0
         RETURN
      ENDIF
C--
C  Set kinematics.
C--
C-- q^2_ZZ.
      Q32MN = ( AM5 + AM6 )**2
      Q32MX = ( RS - ( AM3 + AM4 ) )**2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(8)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
C-- Xi_lo.
      XILO = LOG(Q3*(Q3+AM3+AM4)/S)
      XIHI = LOG(1-2*AM3/RS)
      XI    = XILO + (XIHI-XILO)*Z(5)
      RXI   = EXP(XI)
      WAT   = WAT*(XIHI-XILO)
C--
C      DM1   = (AME*AME/S)*RXI*RXI/(1-RXI)
C      DP1   = (AME*AME/S)
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
      FI1   = Z(9)*x2PI
C-- Phi_2.
      FI2   = MOD(Z(10)*x2PI + FI1,x2PI)
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
      XPHASE(2,3) = -1 + 2*Z(11)
      XPHASE(3,3) = Z(12)*x2PI
C--
      WAT   = WAT*2**(NP-5)*(x2PI)**(NP-3)
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
      XTT1    = 2*SQRT(PV(1,5)**2+PV(2,5)**2)/ROOTS
C--
      XTT2    = 2*SQRT(PV(1,6)**2+PV(2,6)**2)/ROOTS
C--
      XNU     = 2*PV(0,3)/ROOTS
      XNB     = 2*PV(0,4)/ROOTS
      XTP     = 2*PV(0,5)/ROOTS
      XTB     = 2*PV(0,6)/ROOTS
      CSTP    = XPHASE(2,3)
      FITP    = XPHASE(3,3)
      DFI     = MOD(FI1-FI2+x2PI,x2PI)
C--
      CALL XHFILL(  1, DBLE(X1)          ,FUNC )
      CALL XHFILL(  2, DBLE(CS1)         ,FUNC )
      CALL XHFILL(  3, DBLE(FI1)*RDDG    ,FUNC )
      CALL XHFILL(  4, DBLE(X2)          ,FUNC )
      CALL XHFILL(  5, DBLE(CS2)         ,FUNC )
      CALL XHFILL(  6, DBLE(DFI)*RDDG    ,FUNC )
      CALL XHFILL(  7, DBLE(Q3/ROOTS)    ,FUNC )
      CALL XHFILL(  8, DBLE(CSTP)        ,FUNC )
      CALL XHFILL(  9, DBLE(FITP)*RDDG   ,FUNC )
      CALL XHFILL( 10, DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 11, DBLE(XHCMB+.1)    ,FUNC )
      CALL XHFILL( 12, DBLE(ET1)         ,FUNC )
      CALL XHFILL( 13, DBLE(ET2)         ,FUNC )
      CALL XHFILL( 14, DBLE(XTT1)        ,FUNC )
      CALL XHFILL( 15, DBLE(XTT2)        ,FUNC )
      CALL XHFILL( 16, DBLE(XNU)         ,FUNC )
      CALL XHFILL( 17, DBLE(XNB)         ,FUNC )
      CALL XHFILL( 18, DBLE(XTP)         ,FUNC )
      CALL XHFILL( 19, DBLE(XTB)         ,FUNC )
      CALL XHFILL( 20, DBLE(Q3)          ,FUNC )
      CALL XHFILL( 21, DBLE(Q3)          ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
