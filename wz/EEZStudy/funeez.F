C* (Update Record)
C*   95/04/26  K.Fujii	FUNC for e+ e- --> e+ e- Z full amplitudes.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "EEZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6, NHEL = 16 )
      REAL   *8   SG, PV(0:3,NP)
      REAL   *8   EBM, EMBM, EPBM, EMBM0, EPBM0, PEM, PEP
C--
      REAL   *8   QV(0:3,20)
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,-1,   -1, -1,   -1,+1, 
     .             -1,-1,   -1, -1,   +1,-1, 
     .             -1,-1,   -1, +1,   -1,+1, 
     .             -1,-1,   -1, +1,   +1,-1, 
     .             -1,-1,   +1, -1,   -1,+1, 
     .             -1,-1,   +1, -1,   +1,-1, 
     .             -1,-1,   +1, +1,   -1,+1, 
     .             -1,-1,   +1, +1,   +1,-1, 
     .             -1,+1,   -1, -1,   -1,+1, 
     .             -1,+1,   -1, -1,   +1,-1, 
     .             -1,+1,   -1, +1,   -1,+1, 
     .             -1,+1,   -1, +1,   +1,-1, 
     .             -1,+1,   +1, -1,   -1,+1, 
     .             -1,+1,   +1, -1,   +1,-1, 
     .             -1,+1,   +1, +1,   -1,+1, 
     .             -1,+1,   +1, +1,   +1,-1 /
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, -1,   -1,+1, 
     .             +1,-1,   -1, -1,   +1,-1, 
     .             +1,-1,   -1, +1,   -1,+1, 
     .             +1,-1,   -1, +1,   +1,-1, 
     .             +1,-1,   +1, -1,   -1,+1, 
     .             +1,-1,   +1, -1,   +1,-1, 
     .             +1,-1,   +1, +1,   -1,+1, 
     .             +1,-1,   +1, +1,   +1,-1, 
     .             +1,+1,   -1, -1,   -1,+1, 
     .             +1,+1,   -1, -1,   +1,-1, 
     .             +1,+1,   -1, +1,   -1,+1, 
     .             +1,+1,   -1, +1,   +1,-1, 
     .             +1,+1,   +1, -1,   -1,+1, 
     .             +1,+1,   +1, -1,   +1,-1, 
     .             +1,+1,   +1, +1,   -1,+1, 
     .             +1,+1,   +1, +1,   +1,-1 /
C--
      DATA NCALL /  0 /
      EXTERNAL   XHFILL, DHFILL
      DATA EPSL  / 5.D0 /
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
C--
         GMZ   = GMZTOT
         AMGZ  = AMZ*GMZ
         AMZ2  = AMZ*AMZ
         AMGZ2 = AMGZ*AMGZ
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895E0
C--
         XZ1MN = BRVDK(IMDMN-1,2)
         XZ1MX = BRVDK(IMDMX  ,2)
C--
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
C      ( 8) : m(Z)**2
C      ( 9) : phi_2 - phi_1
C      (10) : phi_1
C      (11) : cos_f in Z rest frame.
C      (12) : phi_f in Z rest frame.
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
C--
      IDPTCL(1, 3) = 1
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 1
      IDPTCL(2, 4) = 2
      IDPTCL(3, 4) = 1
C-- Z.
      XZ1 = XZ1MN + Z(13)*(XZ1MX-XZ1MN)
      DO 200 IMD = IMDMN, IMDMX
         IF ( XZ1.LT.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,5) = MDVDK(1,IMD,2)
      IDPTCL(2,5) = MDVDK(2,IMD,2)
      IDPTCL(3,5) = MDVDK(3,IMD,2)
      IDPTCL(1,6) = MDVDK(1,IMD,2)
      IDPTCL(2,6) = MDVDK(2,IMD,2)
      IDPTCL(3,6) = MDVDK(3,IMD,2)
C--
      BRZ1 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      WAT  = WAT*(XZ1MX-XZ1MN)/BRZ1
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
C  Check if RS is sufficient.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
C--
      IF ( RS.LE.AM3+AM4+AM5+AM6+EPSL ) THEN
         FUNC = 0
         RETURN
      ENDIF      
C--
C  Set kinematics.
C--
C-- q^2_Z.
      Q3MN  = AM5 + AM6 + EPSL
      Q3MX  = RS - ( AM3 + AM4 )
      TH3MN = ATAN((Q3MN-AMZ)*(Q3MN+AMZ)/AMGZ)
      TH3MX = ATAN((Q3MX-AMZ)*(Q3MX+AMZ)/AMGZ)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(8)
      Q32   = AMGZ*TAN(TH3) + AMZ2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Xi.
      XILO  = 1/( 1 - 2*AM3/RS )
      XIHI  = S/(Q3*(Q3+AM3+AM4))
      XI    = XILO + (XIHI-XILO)*Z(5)
      RXI   = 1/XI
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
      FI1   = Z(10)*x2PI
      WAT   = WAT*x2PI
C-- Phi_2.
      FI2   = MOD( Z(9)*x2PI+FI1, x2PI )
      WAT   = WAT*x2PI
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
      WAT   = WAT*2*x2PI
C--
C  Calculate differential cross sections.
C--
      CALL SGEEZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
      CALL UVCOPY(8*NP,PVOUT(0,1),PV(0,1))
C--
      EXP2ET  = EXP(2*ET1)
      SH1     = MIN( SQRT( (1+DP1+DM1)/(1+  EXP2ET) - DM1 ), 1. )
      CH1     = MIN( SQRT( (1+DP1+DM1)/(1+1/EXP2ET) - DP1 ), 1. )
      SN1     = 2*SH1*CH1
      CS1     = (1+DP1+DM1)*TANH(ET1) - DP1 + DM1
C--
      EXP2ET  = EXP(2*ET2)
      SH2     = MIN( SQRT( (1+DP2+DM2)/(1+  EXP2ET) - DM2 ), 1. )
      CH2     = MIN( SQRT( (1+DP2+DM2)/(1+1/EXP2ET) - DP2 ), 1. )
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
      CALL PSUMxx(PV(0,5),PV(0,6),QV(0,3))
      XTZ     = 2*SQRT(QV(1,3)**2+QV(2,3)**2)/ROOTS
C--
      XEM     = 2*PV(0,3)/ROOTS
      XEP     = 2*PV(0,4)/ROOTS
      XZ0     = 2*QV(0,3)/ROOTS
      CSFD    = XPHASE(2,3)
      FIFD    = XPHASE(3,3)
      DFI     = MOD(FI2-FI1+x2PI,x2PI)
C--
      CALL XHFILL( 'h01', DBLE(X1)          ,FUNC )
      CALL XHFILL( 'h02', DBLE(CS1)         ,FUNC )
      CALL XHFILL( 'h03', DBLE(FI1)*RDDG    ,FUNC )
      CALL XHFILL( 'h04', DBLE(X2)          ,FUNC )
      CALL XHFILL( 'h05', DBLE(CS2)         ,FUNC )
      CALL XHFILL( 'h06', DBLE(DFI)*RDDG    ,FUNC )
      CALL XHFILL( 'h07', DBLE(Q3)          ,FUNC )
      CALL XHFILL( 'h08', DBLE(CSFD)        ,FUNC )
      CALL XHFILL( 'h09', DBLE(FIFD)*RDDG   ,FUNC )
      CALL XHFILL( 'h10', DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 'h11', DBLE(IMD+  .1D0)  ,FUNC )
      CALL XHFILL( 'h12', DBLE(XHCMB+.1D0)  ,FUNC )
      CALL XHFILL( 'h13', DBLE(ET1)         ,FUNC )
      CALL XHFILL( 'h14', DBLE(ET2)         ,FUNC )
      CALL XHFILL( 'h15', DBLE(XTZ)         ,FUNC )
      CALL XHFILL( 'h16', DBLE(XEM)         ,FUNC )
      CALL XHFILL( 'h17', DBLE(XEP)         ,FUNC )
      CALL XHFILL( 'h18', DBLE(XZ0)         ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
