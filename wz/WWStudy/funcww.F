C* (Update Record)
C*  95/03/12  K.Fujii		Converted to unix. 
C*				To avoid problems related to kinematic
C*				boundaries, we introduced a tiny offset
C*				EPSL. Now use INSMCN instead of CONINT.

      REAL FUNCTION FUNC*8(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
#include "WWCNSF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 6, NHEL = 1 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,5)
      REAL   *8   WAT, WATS
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   -1, +1,   +1,-1,
     .             +1,-1,   -1, +1,   +1,-1/
      DATA NCALL / 0 /
      EXTERNAL   XHFILL, DHFILL
C--
      DATA EPSL / 1. /
C>>>
CCDBG       INTEGER*4  JDP(2,NP)
C>>>
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
C--
         GMW   = GMWTOT
         AMGW  = AMW*GMW
         AMW2  = AMW*AMW
         AMGW2 = AMGW*AMGW
C--
         GMZ   = GMZTOT
         AMGZ  = AMZ*GMZ
         AMZ2  = AMZ*AMZ
         AMGZ2 = AMGZ*AMGZ
C--
         AME   = AMF(1,2,1)
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
C     Z( 1) : e- E_spread
C      ( 2) : e- beamstrahlung
C      ( 3) : e+ E_spread
C      ( 4) : e+ beamstrahlung
C      ( 5) : bremsstrahlng
C      ( 6) : helicity for initial states
C      ( 7) : m(W-)**2
C      ( 8) : m(W+)**2
C      ( 9) : cos(theta_W-)
C      (10) : phi_W-
C      (11) : cos(theta_fd)     in W-   rest frame
C      (12) : phi_fd            in W-   rest frame
C      (13) : cos(theta_fd_bar) in W+   rest frame
C      (14) : phi_fd_bar        in W+   rest frame
C      (15) : final state combination.
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
         CALL PROBEB(REAL(Z(2)),REAL(Z(1)),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(REAL(Z(4)),REAL(Z(3)),EBM,EPBM,EPBM0,PEP)
         RS  = 2*SQRT(EMBM*EPBM)
         IF ( RS.LT.5. ) THEN
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
         ZG  = Z(5)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(3)*1.E6
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
C-- W-.
      XW1 = XW1MN + Z(15)*(XW1MX-XW1MN)
      DO 100 IMD = IMD1MN, IMD1MX
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 110
100   CONTINUE
C--
110   IDPTCL(1,3) = MDVDK(2,IMD,1)
      IDPTCL(2,3) = 2
      IDPTCL(3,3) = MDVDK(3,IMD,1)
      IDPTCL(1,4) = MDVDK(1,IMD,1)
      IDPTCL(2,4) = 1
      IDPTCL(3,4) = MDVDK(3,IMD,1)
C--
      BRW1 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD1 = IMD
      WAT  = WAT*(XW1MX-XW1MN)/BRW1
C
C-- W+.
      XW2 = XW2MN + ((XW1-BRVDK(IMD-1,1))/BRW1)*(XW2MX-XW2MN)
      DO 200 IMD = IMD2MN, IMD2MX
         IF ( XW2.LT.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,5) = MDVDK(2,IMD,1)
      IDPTCL(2,5) = 2
      IDPTCL(3,5) = MDVDK(3,IMD,1)
      IDPTCL(1,6) = MDVDK(1,IMD,1)
      IDPTCL(2,6) = 1
      IDPTCL(3,6) = MDVDK(3,IMD,1)
C--
      BRW2 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD2 = IMD
      WAT  = WAT*(XW2MX-XW2MN)/BRW2
C--
C  Select electron helicity.
C--
      HLM   = (1-POLE)/2
      IF ( Z(6).LT.HLM ) THEN
         IOFF = 0
      ELSE
         IOFF = NHEL
      ENDIF
C--
C  Select final state helicities.
C--
      JHL  = 1
      ICMB = JHL + IOFF
      CALL UVCOPY(NP,IHEL(1,ICMB),IHLCMB)
      WATS  = WAT
C--
C  Set invariant masses squared.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
C-- W-.
      Q2MN  = AM3 + AM4 + EPSL
      Q2MX  = MAX( RS - (AM5+AM6) - EPSL, Q2MN )
      TH2MN = ATAN((Q2MN-AMW)*(Q2MN+AMW)/AMGW)
      TH2MX = ATAN((Q2MX-AMW)*(Q2MX+AMW)/AMGW)
      TH2   = TH2MN + (TH2MX-TH2MN)*Z(7)
      Q22   = AMGW*TAN(TH2) + AMW2
      Q2    = SQRT(Q22)
      WAT   = WAT*(TH2MX-TH2MN)*( (Q22-AMW2)**2 + AMGW2 )/AMGW
C-- W+.
      Q3MN  = AM5 + AM6 + EPSL
      Q3MX  = MAX( RS - Q2 - EPSL, Q3MN )
      TH3MN = ATAN((Q3MN-AMW)*(Q3MN+AMW)/AMGW)
      TH3MX = ATAN((Q3MX-AMW)*(Q3MX+AMW)/AMGW)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(8)
      Q32   = AMGW*TAN(TH3) + AMW2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMW2)**2 + AMGW2 )/AMGW
C--
CZEROWIDTH       Q22   = AMW2
CZEROWIDTH       Q32   = AMW2
CZEROWIDTH       Q2    = AMW
CZEROWIDTH       Q3    = AMW
C--
CZEROWIDTH       WAT   = WATS*(xPI*AMW*GMW)**2
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q22
      XPHASE(1,3) = Q32
C--
C  Angular variables.
C--
      DO 400 IB = 1, 3
         XPHASE(2,IB) = -1 + 2*Z(2*IB+7)
         XPHASE(3,IB) = x2PI  *Z(2*IB+8)
         WAT = WAT*2*x2PI
400   CONTINUE
C>>>
CCDBG       DO 665 K = 1, NP
CCDBG          JDP(1,K) = 100*IDPTCL(1,K)+10*IDPTCL(2,K)+IDPTCL(3,K)
CCDBG          JDP(2,K) = IHLCMB(K)
CCDBG 665   CONTINUE
CCDBG       PRINT *, ' -------    '
CCDBG       PRINT *, ' IDP, IHL = ', ((JDP(K,KK),K=1,2),KK=1,NP)
CCDBG       DO 666 KK = 1, NP-3
CCDBG          PRINT *, ' XPHASE(*,', KK, ') = ', (XPHASE(K,KK),K=1,3)
CCDBG 666   CONTINUE
C>>>
C--
C  Calculate differential cross sections.
C--
      CALL SGWWF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG   = SG*QED
      FUNC = SG*WAT
C>>>
CCDBG       PRINT *, ' SG = ', SG,' FUNC = ', FUNC
C>>>
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
      IF ( ABS(QV(3,1)).GT.1.E-3 ) THEN
C>>>
CCDBG       PRINT *, ' QV = ', (QV(K,1),K=0,3)
C>>>
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
      RDDG  = 180/xPI
      COSW  = XPHASE(2,1)
      PHIW  = XPHASE(3,1)
      COSD1 = XPHASE(2,2)
      PHID1 = XPHASE(3,2)
      COSD2 = XPHASE(2,3)
      PHID2 = XPHASE(3,3)
      CALL XHFILL( 'h01', DBLE(COSW)        ,FUNC )
      CALL XHFILL( 'h02', DBLE(PHIW)*RDDG   ,FUNC )
      CALL XHFILL( 'h03', DBLE(Q2)          ,FUNC )
      CALL XHFILL( 'h04', DBLE(Q3)          ,FUNC )
      CALL XHFILL( 'h05', DBLE(COSD1)       ,FUNC )
      CALL XHFILL( 'h06', DBLE(PHID1)*RDDG  ,FUNC )
      CALL XHFILL( 'h07', DBLE(COSD2)       ,FUNC )
      CALL XHFILL( 'h08', DBLE(PHID2)*RDDG  ,FUNC )
      CALL XHFILL( 'h09', DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 'h10', DBLE(IMD1+.1)     ,FUNC )
      CALL XHFILL( 'h11', DBLE(IMD2+.1)     ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
