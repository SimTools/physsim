C* (Update Record)
C*   2010/10/06  K.Fujii	FUNC for e+ e- --> X0_1 X0_1 full amplitudes.
C*
C

      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "XNXNCN.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6, NHEL = 4 )
      REAL   *8   SG, PV(0:3,NP), AMX(2), GMX(2), AM(NP)
C--
      REAL   *8   AMSF(3,2,2), GTSF(3,2,2)
C--
      REAL   *8   QV(0:3,30)
      REAL   *8   WAT, WT
C--
      INTEGER*4   IHEL(NP,NHEL,2)
      DATA        (( IHEL(I,J,1), I=1,NP ), J=1,NHEL )  /
     .                   -1,+1,   -1, 0,  -1, 0,
     .                   -1,+1,   -1, 0,  +1, 0,
     .                   -1,+1,   +1, 0,  -1, 0,
     .                   -1,+1,   +1, 0,  +1, 0/
      DATA        (( IHEL(I,J,2), I=1,NP ), J=1,NHEL )  /
     .                   +1,-1,   -1, 0,  -1, 0,
     .                   +1,-1,   -1, 0,  +1, 0,
     .                   +1,-1,   +1, 0,  -1, 0,
     .                   +1,-1,   +1, 0,  +1, 0/
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
#if 0
         IX     = 2
         JX     = 2
#else
         IX     = 1
         JX     = 1
#endif
         IF (IX.EQ.JX) THEN
            FACT = 0.5D0
         ELSE
            FACT = 1
         ENDIF
         AMX(1) = SZM (IX)
         GMX(1) = GMSZ(IX)
         AMX(2) = SZM (JX)
         GMX(2) = GMSZ(JX)
C--
         AME    = AMF(1,2,1)
C--
C  3rd generation sfermion masses and widths.
C    (1,2,3,4,5,6,7) = (snuL,seL,seR,suL,suR,sdL,sdR)
C--
C  Ignore D_L and U_2, D_2, assuming that they are heavier
C  than X0_j.
C--
         AMSF(1,1,1) = SFM(1)
         AMSF(1,2,1) = SFM(3)
         AMSF(1,1,2) = SFM(5)
         AMSF(1,2,2) = SFM(7)
C--
         AMSF(2,1,1) = SFM(1)
         AMSF(2,2,1) = SFM(3)
         AMSF(2,1,2) = SFM(5)
         AMSF(2,2,2) = SFM(7)
C--
         AMSF(3,1,1) = SFM3M(1)
         AMSF(3,2,1) = SFM3M(2)
         AMSF(3,1,2) = SFM3M(4)
         AMSF(3,2,2) = SFM3M(6)
C--
         GTSF(1,1,1) = GMSF(1)
         GTSF(1,2,1) = GMSF(3)
         GTSF(1,1,2) = GMSF(5)
         GTSF(1,2,2) = GMSF(7)
C--
         GTSF(2,1,1) = GMSF(1)
         GTSF(2,2,1) = GMSF(3)
         GTSF(2,1,2) = GMSF(5)
         GTSF(2,2,2) = GMSF(7)
C--
         GTSF(3,1,1) = GMSF3(1)
         GTSF(3,2,1) = GMSF3(2)
         GTSF(3,1,2) = GMSF3(4)
         GTSF(3,2,2) = GMSF3(6) 
C--
         xALF0 = 1/137.0359895D0
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
C      ( 4) : XN1 decay mode
C      ( 5) : e- helicity
C      ( 6) : m(X0i)**2
C      ( 7) : m(X0j)**2
C      ( 8) : cos(theta_X0i)
C      ( 9) : phi_X0i
C      (10) : cos(theta_a)     in X0i rest frame
C      (11) : phi_a            in X0i rest frame
C      (12) : cos(theta_a)     in X0j rest frame
C      (13) : phi_a            in X0j rest frame
C      (14) : E_beam spread of E-
C      (15) : E_beam spread of E+
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
         call jsfbeamgend(z(14), z(15), z(1), z(2), ebm, 
     >      embm, epbm, embm0, epbm0, dpdebm)

         RS  = 2*SQRT(EMBM*EPBM)
         IF ( RS.LT.AMX(IX)+AMX(JX) ) THEN
            FUNC = 0
            RETURN
         ENDIF
         BMENGY(0) = EBM
         BMENGY(1) = EMBM0
         BMENGY(2) = EPBM0
         BMENGY(3) = EMBM
         BMENGY(4) = EPBM
cc         DPDEBM    = PEM*PEP
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
         IZE = Z(9)*1.D8
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
C--
C  Set initial and final state particles:
C     IDPTCL(1,IP) = generation number (1~3)
C           (2,IP) = isospin = (1,2) = (up,down)
C           (3,IP) = (1,2) = (lepton,quark)
C     with
C        IP = 1: e-
C           = 2: e+
C           = 3: tau-/+
C           = 4: stau+/-
C           = 5: tau+/-
C           = 6: stau-/+
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
C-- 1st X0.
C     IMD=(1,2)=((tau-,stau+),(tau+,stau-))
C
      XZ1 = Z(4)
      IF ( IMD1MX.LE.1 ) THEN
         IMD1 = 1
         XZ2  = XZ1
      ELSE IF ( IMD1MN.GE.2 ) THEN
         IMD1 = 2
         XZ2  = XZ1
      ELSE
         IF ( XZ1.LT.0.5D0 ) THEN
            IMD1 = 1
            XZ2  = XZ1*2
         ELSE
            IMD1 = 2
            XZ2  = (XZ1-0.5D0)*2
         ENDIF
         WAT  = WAT*2
      ENDIF
C--
      IDPTCL(1, 3) = 3
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 3
      IDPTCL(2, 4) = 2*(-1)**(IMD1-1)
      IDPTCL(3, 4) = 1
C--
C-- 2nd X0.
      IF ( IMD2MX.LE.1 ) THEN
         IMD2 = 1
      ELSE IF ( IMD2MN.GE.2 ) THEN
         IMD2 = 2
      ELSE
         IF ( XZ2.LT.0.5D0 ) THEN
            IMD2 = 1
         ELSE
            IMD2 = 2
         ENDIF
         WAT  = WAT*2
      ENDIF
C--
      IDPTCL(1, 5) = 3
      IDPTCL(2, 5) = 2
      IDPTCL(3, 5) = 1
      IDPTCL(1, 6) = 3
      IDPTCL(2, 6) = 2*(-1)**(IMD2-1)
      IDPTCL(3, 6) = 1
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( Z(5).LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*Z(5)/HLM + 1
      ELSE
         ICMB = 2
         JCMB = NHEL*(Z(5)-HLM)/(1-HLM) + 1
      ENDIF
      WAT = WAT*NHEL
      JCMB = MIN(JCMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
C--
C  Select kinematic variables.
C     (1,2,3,4,5,6) = (e-,e+, tau-/+,stau+/-, tau+/-,stau-/+)
C--
      AM1   = AMF (IDPTCL(1,1),IDPTCL(2,1),IDPTCL(3,1))
      AM2   = AMF (IDPTCL(1,2),IDPTCL(2,2),IDPTCL(3,2))
      AM3   = AMF (IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMSF(IDPTCL(1,4),IABS(IDPTCL(2,4)),IDPTCL(3,4))
      AM5   = AMF (IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMSF(IDPTCL(1,6),IABS(IDPTCL(2,6)),IDPTCL(3,6))
      AM(1) = AM1
      AM(2) = AM2
      AM(3) = AM3
      AM(4) = AM4
      AM(5) = AM5
      AM(6) = AM6
C--
C  Set phase space variables.
C-- 
C-- 1st X0.
      Q3MN   = AM3 + AM4
      Q3MX   = RS - (AM5+AM6)
      AMXM   = AMX(IX)
      AMXM2  = AMXM*AMXM
      AMGXM  = AMXM*GMX(IX)
      AMGXM2 = AMGXM*AMGXM
      TH3MN  = ATAN((Q3MN-AMXM)*(Q3MN+AMXM)/AMGXM)
      TH3MX  = ATAN((Q3MX-AMXM)*(Q3MX+AMXM)/AMGXM)
      TH3    = TH3MN + (TH3MX-TH3MN)*Z(6)
      Q32    = AMGXM*TAN(TH3) + AMXM2
      Q3     = SQRT(Q32)
      WAT    = WAT*(TH3MX-TH3MN)*( (Q32-AMXM2)**2 + AMGXM2 )/AMGXM
C-- 2nd X0.
      Q4MN   = AM5 + AM6
      Q4MX   = RS - Q3
      AMXP   = AMX(JX)
      AMXP2  = AMXP*AMXP
      AMGXP  = AMXP*GMX(JX)
      AMGXP2 = AMGXP*AMGXP
      TH4MN  = ATAN((Q4MN-AMXP)*(Q4MN+AMXP)/AMGXP)
      TH4MX  = ATAN((Q4MX-AMXP)*(Q4MX+AMXP)/AMGXP)
      TH4    = TH4MN + (TH4MX-TH4MN)*Z(7)
      Q42    = AMGXP*TAN(TH4) + AMXP2
      Q4     = SQRT(Q42)
      WAT    = WAT*(TH4MX-TH4MN)*( (Q42-AMXP2)**2 + AMGXP2 )/AMGXP
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q32
      XPHASE(1,3) = Q42
C--
C  Angular variables.
C--
      DO 300 IB = 1, 3
         XPHASE(2,IB) = -1 + 2*Z(2*IB+6)
         XPHASE(3,IB) = x2PI  *Z(2*IB+7)
         WAT = WAT*2*x2PI
300   CONTINUE

C--
C  Calculate differential cross sections.
C--
      CALL SGNXXF(IDPTCL,IHLCMB,AM,XPHASE,SG,PV)
      SG     = SG*QED*FACT
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
      CALL UVCOPY(2*4*NP,PVOUT(0,1),PV(0,1))
C--
      EBM   = RS/2
      S     = RS*RS
C--      
      COST   = XPHASE(2,1)
      PHIT   = XPHASE(3,1)
      COSBB1 = XPHASE(2,2)
      PHIBB1 = XPHASE(3,2)
      COSBB2 = XPHASE(2,3)
      PHIBB2 = XPHASE(3,3)
      COMB   = JCMB + NHEL*(ICMB-1) + .1D0

      QXI2   = XPHASE(1,2)
      QXJ2   = XPHASE(1,3)
C--      
      CALL XHFILL( 'h01', COST               ,FUNC )
      CALL XHFILL( 'h02', PHIT*xRD2DG        ,FUNC )
      CALL XHFILL( 'h03', COSBB1             ,FUNC )
      CALL XHFILL( 'h04', PHIBB1*xRD2DG      ,FUNC )
      CALL XHFILL( 'h05', SQRT(QXI2)         ,FUNC )
      CALL XHFILL( 'h06', COSBB2             ,FUNC )
      CALL XHFILL( 'h07', PHIBB2*xRD2DG      ,FUNC )
      CALL XHFILL( 'h08', SQRT(QXJ2)         ,FUNC )
      CALL XHFILL( 'h09', RS/ROOTS           ,FUNC )
      CALL XHFILL( 'h10', COMB               ,FUNC )
      CALL XHFILL( 'h11', IMD1+.1D0          ,FUNC )
      CALL XHFILL( 'h11', IMD2+2.1D0         ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
