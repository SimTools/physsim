C* (Update Record)
C*   95/04/26  K.Fujii	FUNC for e+ e- --> nu nubar Z full amplitudes.
C*
      REAL FUNCTION FUNC*8(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "NNZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6, NHEL = 2 )
      REAL   *8   SG, PV(0:3,NP)
C--
      REAL   *8   QV(0:3,20)
      REAL   *8   Q3, Q32, Q3MN, Q3MX, TH3MN, TH3MX, TH3
      REAL   *8   WAT
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,+1,   -1, +1,   +1,-1,
     .             -1,+1,   -1, +1,   -1,+1/
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, +1,   +1,-1,
     .             +1,-1,   -1, +1,   -1,+1/
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
         AMW2  = AMW*AMW
         GMZ   = GMZTOT
         AMGZ  = AMZ*GMZ
         AMZ2  = AMZ*AMZ
         AMGZ2 = AMGZ*AMGZ
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895D0
C--
         XZ1MN = BRVDK(IMDMN-1,2)
         XZ1MX = BRVDK(IMDMX  ,2)
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
C      ( 5) : xi   (E_nu)
C      ( 6) : m(Z)**2
C      ( 7) : eta  (cos_nu)
C      ( 8) : zeta (cos_nb) in (nb Z) frame
C      ( 9) : phi_nu
C      (10) : phi_nb  	   in (nb Z) frame
C      (11) : cos_fb       in Z      frame
C      (12) : phi_f        in Z      frame
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
      IDPTCL(2, 3) = 1
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 1
      IDPTCL(2, 4) = 1
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
C  Set kinematics.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
C-- q^2_Z.
      Q3MN  = AM5 + AM6
      Q3MX  = RS - ( AM3 + AM4 )
      TH3MN = ATAN((Q3MN-AMZ)*(Q3MN+AMZ)/AMGZ)
      TH3MX = ATAN((Q3MX-AMZ)*(Q3MX+AMZ)/AMGZ)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(6)
      Q32   = AMGZ*TAN(TH3) + AMZ2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Xi_lo.
      XILO  = LOG( (Q3-AM3)*(Q3+AM3)/S )
      XIHI  = LOG( 1 - 2*AM3/RS )
      XI    = XILO + (XIHI-XILO)*Z(5)
      WAT   = WAT*(XIHI-XILO)
      EXPXI = EXP(XI)
C-- Eta.
      DM1   = AMW2/S
      DP1   = 1
      ETLO  = -LOG( (1+DM1)/DP1 )/2
      ETHI  =  LOG( (1+DP1)/DM1 )/2
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
      CALL SGNNZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
C  Fill hists. and plots.
C--
      RDDG   = 180/xPI
      QWN2   = S*EXPXI + AM3*AM3
      QWN    = SQRT(QWN2)
      EXP2ET = EXP(2*XPHASE(2,1))
      SH1    = MIN( SQRT( (1+DP1+DM1)/(1+EXP2ET) - DM1 ), 1.D0 )
      DLT    = 2*AMZ*AMZ/((QWN-AMZ)*(QWN+AMZ))
      OPDZT  = ( 1 + 1/DLT )**XPHASE(2,2)
      SH2    = MIN( SQRT( DLT*( OPDZT - 1 ) ), 1.D0 )
C--
      XE    = 1 + AM3*AM3/S - QWN2/S
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
      XTW     = 2*SQRT(QV(1,5)**2+QV(2,5)**2+QV(3,5)**2)/ROOTS
C--
      XEL     = 2*PV(0,3)/ROOTS
      XNU     = 2*PV(0,4)/ROOTS
      XWP     = 2*QV(0,5)/ROOTS
C--
      CALL XHFILL( 'h01', COSE              ,FUNC )
      CALL XHFILL( 'h02', PHIE*RDDG         ,FUNC )
      CALL XHFILL( 'h03', XE                ,FUNC )
      CALL XHFILL( 'h04', QW                ,FUNC )
      CALL XHFILL( 'h05', COSN              ,FUNC )
      CALL XHFILL( 'h06', PHIN*RDDG         ,FUNC )
      CALL XHFILL( 'h07', CSDA              ,FUNC )
      CALL XHFILL( 'h08', PHID*RDDG         ,FUNC )
      CALL XHFILL( 'h09', RS/ROOTS          ,FUNC )
      CALL XHFILL( 'h10', IMD+.1D0          ,FUNC )
      CALL XHFILL( 'h11', XHCMB+.1D0        ,FUNC )
      CALL XHFILL( 'h12', ETA               ,FUNC )
      CALL XHFILL( 'h13', XTW               ,FUNC )
      CALL XHFILL( 'h14', XEL               ,FUNC )
      CALL XHFILL( 'h15', XNU               ,FUNC )
      CALL XHFILL( 'h16', XWP               ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
