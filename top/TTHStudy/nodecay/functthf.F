C* 95/06/08  K.Fujii    FUNC for e+ e- --> t t_bar H w/o decays.
C*

      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
#include "TTZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 5, NHEL = 4 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,5)
C--
      REAL   *8   WAT, WATS
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   -1, -1, 0,
     .             -1,+1,   -1, +1, 0,
     .             -1,+1,   +1, -1, 0,
     .             -1,+1,   +1, +1, 0,
     .             +1,-1,   -1, -1, 0,
     .             +1,-1,   -1, +1, 0,
     .             +1,-1,   +1, -1, 0,
     .             +1,-1,   +1, +1, 0/
C>>>
CCDBG       INTEGER*4  JDP(2,NP)
C>>>
      DATA NCALL / 0 /
C--
      DATA EPSL  / 0.0 /
C
C========< Entry Point >================================================
C
C>>>
CCDBG        print *, '>>>> Now enters FUNC.'
C>>>
C--
C  Set some variables.
C     ENGM  = variable mass range in units of width.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
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
C      ( 3) : ISR
C      ( 4) : e- helicity
C      ( 5) : helicity combination for final states.
C      ( 6) : m(t-t_bar)**2
C      ( 7) : cos(theta_H)
C      ( 8) : phi_H
C      ( 9) : cos(theta_t)     in t-t_bar rest frame
C      (10) : phi_t            in t-t_bar rest frame
C      (11) : e- beam gaussian spread
C      (12) : e+ beam gaussian spread
C--
      WAT = 1
C--
      QED = 1
C--
C  First decide e- and e+ energy.
C--
      EBM       = ROOTS/2
      RS        = ROOTS
      EMBM      = EBM
      EPBM      = EBM
C--
      BMENGY(0) = EBM
      BMENGY(1) = BMENGY(0)
      BMENGY(2) = BMENGY(0)
      BMENGY(3) = BMENGY(0)
      BMENGY(4) = BMENGY(0)
      DPDEBM    = 1
C--      
      IF ( ISRBM.EQ.3 ) THEN
         CALL PROBEB(REAL(Z(1)),REAL(Z(11)),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(REAL(Z(2)),REAL(Z(12)),EBM,EPBM,EPBM0,PEP)
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
C  Then decide reduced sqrt(s) after ISR.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(8)*1.E6
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
C--      
      S      = RS*RS
C--
C  Now set initial states.
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
#if 1
C tt
      IDPTCL(1, 3) = 3
      IDPTCL(2, 3) = 1
      IDPTCL(3, 3) = 2
      IDPTCL(1, 4) = 3
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = 2
#else
C bb
      IDPTCL(1, 3) = 3
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 2
      IDPTCL(1, 4) = 3
      IDPTCL(2, 4) = 2
      IDPTCL(3, 4) = 2
#endif
C--
C  Select electron helicity.
C--
      HLM   = (1-POLE)/2
      IF ( Z(4).LT.HLM ) THEN
         IOFF = 0
      ELSE
         IOFF = NHEL
      ENDIF
C--
C  Select final state helicities.
C--
      JHL  = NHEL*Z(5)+1
      JHL  = MIN(JHL,NHEL)
      ICMB = JHL + IOFF
      WAT  = WAT*NHEL
      CALL UVCOPY(NP,IHEL(1,ICMB),IHLCMB)
      WATS  = WAT
C>>>
CCDBG        print *, 'Now calculate kinematics'
C>>>
C--
C  Set .
C--
      AM3   = AMF(IDPTCL(1, 3),IDPTCL(2, 3),IDPTCL(3, 3))
      AM4   = AMF(IDPTCL(1, 4),IDPTCL(2, 4),IDPTCL(3, 4))
      AM5   = AMH
C--
C  Check if RS exceeds the sum of the masses.
C--
      AMSUM = AM3 + AM4 + AM5
      IF ( RS.LE.AMSUM ) THEN
         FUNC = 0
         RETURN
      ENDIF
C--
C  Set invariant masses squared.
C--
C>>>
CCDBG        print *, '--ttbar'
C>>>
C-- t-t_bar.
      Q62MN = (AM3+AM4)**2
      Q62MX = (RS-AM5)**2
      Q62   = Q62MN + (Q62MX-Q62MN)*Z(6)
      WAT   = WAT*(Q62MX-Q62MN)
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q62
C--
C  Angular variables.
C--
      DO 400 IB = 1, 2
         XPHASE(2,IB) = -1 + 2*Z(2*IB+5)
         XPHASE(3,IB) = x2PI  *Z(2*IB+6)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Calculate differential cross sections.
C--
C>>>
CCDBG        print *, 'Going to call SGTTHF'
C>>>
      CALL SGTTHF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG   = SG*QED
      FUNC = SG*WAT
C>>>
CCDBG       DO 665 K = 1, NP
CCDBG          JDP(1,K) = 100*IDPTCL(1,K)+10*IDPTCL(2,K)+IDPTCL(3,K)
CCDBG          JDP(2,K) = IHLCMB(K)
CCDBG 665   CONTINUE
CCDBG       PRINT *, ' -------    '
CCDBG       PRINT *, ' IDP, IHL = ', ((JDP(K,KK),K=1,2),KK=1,NP)
CCDBG       DO 666 KK = 1, 2
CCDBG          PRINT *, ' XPHASE(*,', KK, ') = ', (XPHASE(K,KK),K=1,3)
CCDBG 666   CONTINUE
CCDBG       PRINT *, ' SG = ', SG,' FUNC = ', FUNC
C>>>
      IF ( SG.LE.0. ) RETURN
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
C>>>
CCDBG        print *, 'Hists and plots.'
C>>>
C--
      AMSTT2 = XPHASE(1,2)
      AMTT  = SQRT(AMSTT2)
      RDDG  = 180/xPI
      COSZ  = XPHASE(2,1)
      PHIZ  = XPHASE(3,1)
      COST  = XPHASE(2,2)
      PHIT  = XPHASE(3,2)
      CALL XHFILL(  1, DBLE(COSZ)        ,FUNC )
      CALL XHFILL(  2, DBLE(PHIZ)*RDDG   ,FUNC )
      CALL XHFILL(  3, DBLE(AMTT/ROOTS)  ,FUNC )
      CALL XHFILL(  4, DBLE(COST)        ,FUNC )
      CALL XHFILL(  5, DBLE(PHIT)*RDDG   ,FUNC )
      CALL XHFILL(  6, DBLE(CST)         ,FUNC )
      CALL XHFILL(  7, DBLE(ICMB)        ,FUNC )
C>>>
CCDBG        print *, 'That''s it.'
C>>>
C--
C  That's it.
C--
      RETURN
      END
