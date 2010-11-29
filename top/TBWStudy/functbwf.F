C* 2010/11/29  K.Fujii	Derived from TTStudy/functtf.F
C*
#if 0
#define __ZEROWIDTH__
#endif
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
      INCLUDE 'TBWCNF.inc'
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 8, NHEL = 1 )
      REAL   *8   SG, PV(0:3,NP), QV(0:3,10)
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   +1, +1,-1,   -1, -1,+1,
     .             +1,-1,   +1, +1,-1,   -1, -1,+1/
C>>>
CCDBG       INTEGER*4  JDP(2,NP)
CCDBG       REAL   *8  PVTEST(0:3,10)
C>>>
      DATA EPS   / 0.0D0 /
      DATA NCALL / 0 /
      EXTERNAL   XHFILL, DHFILL
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT( MAX( 1 - 2*(X1+X2) + (X1-X2)**2, 0.D0 ) )
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         GMW   = GMWTOT
         AMGW  = AMW*GMW
         AMW2  = AMW*AMW
         AMGW2 = AMGW*AMGW
C--
         AME   = AMF(1,2,1)
         xALF0 = 1/137.0359895D0
C--
         XW1MN = BRVDK(IMD1MN-1,1)
         XW1MX = BRVDK(IMD1MX  ,1)
         XW2MN = BRVDK(IMD2MN-1,1)
         XW2MX = BRVDK(IMD2MX  ,1)
      ENDIF
C--
      AMT   = AMF(3,1,2)
      GMT   = GMF(3,1,2)
      AMGT  = AMT*GMT
      AMT2  = AMT*AMT
      AMGT2 = AMGT*AMGT
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
C      ( 4) : e- helicity
C      ( 5) : m(t_bar)**2
C      ( 6) : m(t)**2
C      ( 7) : m(W-)**2
C      ( 8) : m(W+)**2
C      ( 9) : cos(theta_t)
C      (10) : phi_t
C      (11) : cos(theta_q_bar) in t_bar rest frame
C      (12) : phi_q_bar        in t_bar rest frame
C      (13) : cos(theta_f)     in W- rest frame
C      (14) : phi_f            in W- rest frame
C      (15) : cos(theta_q)     in t rest frame
C      (16) : phi_q            in t rest frame
C      (17) : cos(theta_f_bar) in W+ rest frame
C      (18) : phi_f_bar        in W+ rest frame
C      (19) : final state combination.
C      (20) : e- beam gaussian spread
C      (21) : e+ beam gaussian spread
C--
      WAT = 1
C--
C  First decide e- and e+ energy.
C--
      EBM = ROOTS/2
      RS        = ROOTS
      EMBM      = EBM
      EPBM      = EBM
      QED       = 1
C--
      BMENGY(0) = EBM
      BMENGY(1) = BMENGY(0)
      BMENGY(2) = BMENGY(0)
      BMENGY(3) = BMENGY(0)
      BMENGY(4) = BMENGY(0)
      DPDEBM    = 1
C--      
      IF ( ISRBM.EQ.3 ) THEN
         call jsfbeamgenr(Z(1), Z(2), Z(20), Z(21), ebm,
     >                    embm, epbm, embm0, epbm0, dpdebm )
         RS  = 2*SQRT(EMBM*EPBM)
         IF ( RS.LT.5.D0 ) THEN
            FUNC = 0
            RETURN
         ENDIF
         BMENGY(0) = EBM
         BMENGY(1) = EMBM0
         BMENGY(2) = EPBM0
         BMENGY(3) = EMBM
         BMENGY(4) = EPBM
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
         IZE = Z(19)*1.D6
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
C--      
      S      = RS*RS
C--
C  Set initial states.
C--
      IDPTCL(1, 1) = 1
      IDPTCL(2, 1) = 2
      IDPTCL(3, 1) = 1
      IDPTCL(1, 2) = 1
      IDPTCL(2, 2) = 2
      IDPTCL(3, 2) = 1
C--
C  Then select final states.
C     t --> b + W only.
C--
      IDPTCL(1, 3) = 3
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 2
C>>>
CTBWOTSW
CTBWOTSW       IDPTCL(1, 3) = 2
CTBWOTSW       IDPTCL(2, 3) = 2
CTBWOTSW       IDPTCL(3, 3) = 2
C>>>
      IDPTCL(1, 6) = 3
      IDPTCL(2, 6) = 2
      IDPTCL(3, 6) = 2
C
C-- W from t-bar.
      XW1 = XW1MN + Z(19)*(XW1MX-XW1MN)
      DO 200 IMD = IMD1MN, IMD1MX
         IF ( XW1.LE.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1, 4) = MDVDK(1,IMD,1)
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = MDVDK(3,IMD,1)
      IDPTCL(1, 5) = MDVDK(2,IMD,1)
      IDPTCL(2, 5) = 2
      IDPTCL(3, 5) = MDVDK(3,IMD,1)
C--
      BRW1 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD1 = IMD
      WAT  = WAT*(XW1MX-XW1MN)/BRW1
C
C-- W from t.
      XW2 = XW2MN + ((XW1-BRVDK(IMD-1,1))/BRW1)*(XW2MX-XW2MN)
      DO 300 IMD = IMD2MN, IMD2MX
         IF ( XW2.LE.BRVDK(IMD,1) )              GO TO 310
300   CONTINUE
C--
310   IDPTCL(1, 7) = MDVDK(1,IMD,1)
      IDPTCL(2, 7) = 1
      IDPTCL(3, 7) = MDVDK(3,IMD,1)
      IDPTCL(1, 8) = MDVDK(2,IMD,1)
      IDPTCL(2, 8) = 2
      IDPTCL(3, 8) = MDVDK(3,IMD,1)
C--
      BRW2 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD2 = IMD
      WAT  = WAT*(XW2MX-XW2MN)/BRW2
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( Z(4).LT.HLM ) THEN
         ICMB = 1
      ELSE
         ICMB = 2
      ENDIF
      CALL UVCOPY(NP,IHEL(1,ICMB),IHLCMB(1))
C--
C  Set invariant masses squared.
C--
#ifdef __ZEROWIDTH__
      Q32   = AMT2
      Q42   = AMT2
      Q72   = AMW2
      Q82   = AMW2
      WAT   = WAT*(xPI*AMT*GMT)**2
     .           *(xPI*AMW*GMW)**2
                                        GO TO 350
#endif
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AM7   = AMF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AM8   = AMF(IDPTCL(1,8),IDPTCL(2,8),IDPTCL(3,8))
C-- t_bar.
      Q3MN  = AM3 + AM4 + AM5 + EPS
      Q3MX  = RS - (AM6+AM7+AM8+EPS)
#if 0
      TH3MN = ATAN((Q3MN-AMT)*(Q3MN+AMT)/AMGT)
      TH3MX = ATAN((Q3MX-AMT)*(Q3MX+AMT)/AMGT)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(5)
      Q32   = AMGT*TAN(TH3) + AMT2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMT2)**2 + AMGT2 )/AMGT
#else
      Q32MN = Q3MN*Q3MN
      Q32MX = Q3MX*Q3MX
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(5)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
#endif
C-- t.
      Q4MN  = AM6 + AM7 + AM8 + EPS
      Q4MX  = RS - Q3
#if 0
      TH4MN = ATAN((Q4MN-AMT)*(Q4MN+AMT)/AMGT)
      TH4MX = ATAN((Q4MX-AMT)*(Q4MX+AMT)/AMGT)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(6)
      Q42   = AMGT*TAN(TH4) + AMT2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMT2)**2 + AMGT2 )/AMGT
#else
      Q42MN = Q4MN*Q4MN
      Q42MX = Q4MX*Q4MX
      Q42   = Q42MN + (Q42MX-Q42MN)*Z(6)
      Q4    = SQRT(Q42)
      WAT   = WAT*(Q42MX-Q42MN)
#endif
C-- W-.
      Q7MN  = AM4 + AM5
      Q7MX  = Q3 - AM3
      TH7MN = ATAN((Q7MN-AMW)*(Q7MN+AMW)/AMGW)
      TH7MX = ATAN((Q7MX-AMW)*(Q7MX+AMW)/AMGW)
      TH7   = TH7MN + (TH7MX-TH7MN)*Z(7)
      Q72   = AMGW*TAN(TH7) + AMW2
      WAT   = WAT*(TH7MX-TH7MN)*( (Q72-AMW2)**2 + AMGW2 )/AMGW
C-- W+.
      Q8MN  = AM7 + AM8
      Q8MX  = Q4 - AM6
      TH8MN = ATAN((Q8MN-AMW)*(Q8MN+AMW)/AMGW)
      TH8MX = ATAN((Q8MX-AMW)*(Q8MX+AMW)/AMGW)
      TH8   = TH8MN + (TH8MX-TH8MN)*Z(8)
      Q82   = AMGW*TAN(TH8) + AMW2
      WAT   = WAT*(TH8MX-TH8MN)*( (Q82-AMW2)**2 + AMGW2 )/AMGW
C--
350   XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q32
      XPHASE(1,3) = Q72
      XPHASE(1,4) = Q42
      XPHASE(1,5) = Q82
C--
C  Angular variables.
C--
      DO 400 IB = 1, 5
         XPHASE(2,IB) = -1 + 2*Z(2*IB+7)
         XPHASE(3,IB) = x2PI  *Z(2*IB+8)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Calculate differential cross sections.
C--
      CALL SGTBWF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG     = SG*QED
      FUNC   = WAT*SG
      DSGMDX = SG
C>>>
CCDBG       DO 665 K = 1, NP
CCDBG          JDP(1,K) = 100*IDPTCL(1,K)+10*IDPTCL(2,K)+IDPTCL(3,K)
CCDBG          JDP(2,K) = IHLCMB(K)
CCDBG 665   CONTINUE
CCDBG       PRINT *, ' -------    '
CCDBG       PRINT *, ' IDP, IHL = ', ((JDP(K,KK),K=1,2),KK=1,NP)
CCDBG       DO 666 KK = 1, 5
CCDBG          PRINT *, ' XPHASE(*,', KK, ') = ', (XPHASE(K,KK),K=1,3)
CCDBG 666   CONTINUE
CCDBG       CALL UVZERO(2*4,PVTEST)
CCDBG       DO 667 KK = 1, 8
CCDBG          PRINT *, ' PV(*,', KK, ') = ', (PV(K,KK),K=0,3)
CCDBG          IF ( KK.GE.3 )
CCDBG      .         CALL PSUMxx(PV(0,KK),PVTEST(0,1),PVTEST(0,1))
CCDBG 667   CONTINUE
CCDBG       PRINT *, ' PVSUM = ', (PVTEST(K,1),K=0,3)
CCDBG       CALL PSUMxx(PV(0,4),PV(0,5),PVTEST(0,7))
CCDBG       CALL PSUMxx(PVTEST(0,7),PV(0,3),PVTEST(0,3))
CCDBG       CALL PSUMxx(PV(0,7),PV(0,8),PVTEST(0,8))
CCDBG       CALL PSUMxx(PVTEST(0,8),PV(0,6),PVTEST(0,4))
CCDBG       CALL RDOTxx(PVTEST(0,7),PVTEST(0,7),AM)
CCDBG       AM = SIGN(SQRT(ABS(AM)),AM)
CCDBG       PRINT *, ' AM, PV(W-) = ', AM, (PVTEST(K,7),K=0,3)
CCDBG       CALL RDOTxx(PVTEST(0,3),PVTEST(0,3),AM)
CCDBG       AM = SIGN(SQRT(ABS(AM)),AM)
CCDBG       PRINT *, ' AM, PV(tB) = ', AM, (PVTEST(K,3),K=0,3)
CCDBG       CALL RDOTxx(PVTEST(0,8),PVTEST(0,8),AM)
CCDBG       AM = SIGN(SQRT(ABS(AM)),AM)
CCDBG       PRINT *, ' AM, PV(W+) = ', AM, (PVTEST(K,8),K=0,3)
CCDBG       CALL RDOTxx(PVTEST(0,4),PVTEST(0,4),AM)
CCDBG       AM = SIGN(SQRT(ABS(AM)),AM)
CCDBG       PRINT *, ' AM, PV(t ) = ', AM, (PVTEST(K,4),K=0,3)
CCDBG       PRINT *, ' SG = ', SG,' FUNC = ', FUNC
C>>>
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
      IF ( ABS(QV(3,1)).GT.1.D-3 ) THEN
C>>>
CCDBG       PRINT *, ' QV = ', (QV(K,1),K=0,3)
C>>>
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(2*4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
      COST  = XPHASE(2,1)
      PHIT  = XPHASE(3,1)
      COSBB = XPHASE(2,2)
      PHIBB = XPHASE(3,2)
      COSB  = XPHASE(2,4)
      PHIB  = XPHASE(3,4)
      PTOP  = (RS/2)*BETA(Q42/RS**2,Q32/RS**2)
      CALL PSUMxx(PVOUT(0,4),PVOUT(0, 5),QV(0,3))
      CALL PSUMxx(QV(0,3)   ,PVOUT(0, 3),QV(0,3))
      CALL PSUMxx(PVOUT(0,7),PVOUT(0, 8),QV(0,4))
      CALL PSUMxx(QV(0,4)   ,PVOUT(0, 6),QV(0,4))
      CALL PSUMxx(QV(0,4)   ,PVOUT(0, 3),QV(0,9))
      CALL RDOTxx(QV(0,9),QV(0,9),QTBB2)
      QTBB = SQRT(QTBB2)

      PT1  = SQRT(QV(1,3)**2+QV(2,3)**2+QV(3,3)**2)
      PT2  = SQRT(QV(1,4)**2+QV(2,4)**2+QV(3,4)**2)
C--
      CALL XHFILL( 'h01', COST              ,FUNC )
      CALL XHFILL( 'h02', PHIT*xRD2DG       ,FUNC )
      CALL XHFILL( 'h03', COSBB             ,FUNC )
      CALL XHFILL( 'h04', PHIBB*xRD2DG      ,FUNC )
      CALL XHFILL( 'h05', COSB              ,FUNC )
      CALL XHFILL( 'h06', PHIB *xRD2DG      ,FUNC )
      CALL XHFILL( 'h07', RS/ROOTS          ,FUNC )
      CALL XHFILL( 'h08', QTBB              ,FUNC )
      CALL XHFILL( 'h09', Q4                ,FUNC )
      CALL XHFILL( 'h10', Q3                ,FUNC )
      CALL XHFILL( 'h11', PTOP              ,FUNC )
      CALL XHFILL( 'h12', PT1               ,FUNC )
      CALL XHFILL( 'h12', PT2               ,FUNC )
      CALL XHFILL( 'h13', IMD1+.1D0         ,FUNC )
      CALL XHFILL( 'h14', IMD2+.1D0         ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
