C*  7/02/90  K.Fujii    FUNC for e+ e- --> t t_bar full amplitudes.
C* 10/13/90  K.Fujii    Takes into account beam and bremsstrahlng
C*                      as well as resonance shape.
C* 10/27/92  K.Fujii    QED correction ( 1+3*BTE/4 ) and
C*                      ( 1 + 2*ALF0/PI * ( PI**2/6 - 1/4 ) )
C*                      included.
C* 95/06/07  K.Fujii	Now use INSMCN instead of CONINT.
C*			Use ISRBM switch.
C*
#if 0
#define __ZEROWIDTH__
#endif
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
      REAL*8   EPS
      DATA     EPS/4.D-2/
C--
C  Generator COMMON's.
C--
#include "TTZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
#include "USRPRM.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 8, NHEL = 1 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,10)
#ifdef WITH_DBD_STANDARD
      REAL*8 PISR(0:3,2)
      REAL*8 WAT
#endif
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   +1, +1,-1,   -1, -1,+1,
     .             +1,-1,   +1, +1,-1,   -1, -1,+1/
C>>>
CCDBG       INTEGER*4  JDP(2,NP)
CCDBG       REAL   *4  PVTEST(0:3,10)
C>>>
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
         xALF0 = 1/137.0359895E0
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
C      ( 3) : ISR for e-
C      ( 4) : ISR for e+ ( or e-/e+ selection when ISR_LLR_Order=0)
C      ( 5) : e- helicity
C      ( 6) : m(t_bar)**2
C      ( 7) : m(t)**2
C      ( 8) : m(W-)**2
C      ( 9) : m(W+)**2
C      (10) : cos(theta_t)
C      (11) : phi_t
C      (12) : cos(theta_q_bar) in t_bar rest frame
C      (13) : phi_q_bar        in t_bar rest frame
C      (14) : cos(theta_f)     in W- rest frame
C      (15) : phi_f            in W- rest frame
C      (16) : cos(theta_q)     in t rest frame
C      (17) : phi_q            in t rest frame
C      (18) : cos(theta_f_bar) in W+ rest frame
C      (19) : phi_f_bar        in W+ rest frame
C      (20) : final state combination.
C      (21) : e- beam gaussian spread
C      (22) : e+ beam gaussian spread
C      (23) : random variable for Pt of ISR photon from e- beam
C      (24) : random variable for Pt of ISR photon from e+ beam
C      (25) : random variable for azimuthal angle of ISR photon from e- beam
C      (26) : random variable for azimuthal angle of ISR photon from e+ beam

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
#ifndef WITH_DBD_STANDARD
      IF ( ISRBM.EQ.3 ) THEN
         call jsfbeamgenr(Z(1), Z(2), Z(20), Z(21), ebm,
     >                    embm, epbm, embm0, epbm0, dpdebm )
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
         wat = wat*dpdebm
      ENDIF
#else
C -- when WITH_DBD_STANDARD is set
      IF ( ISRBM.GE.3 ) THEN
        IF ( ISRBM.EQ.3 ) THEN
          CALL JSFBEAMGENR(Z(1), Z(2), Z(21), Z(22), EBM,
     >                    EMBM, EPBM, EMBM0, EPBM0, DPDEBM )
        ELSE
          CALL WHIZARD_SPECTRUM( ISRBM, ROOTS, Z(1), Z(2),
     >       EBM, DPDEBM, EMBM, EPBM )
          EMBM0=EMBM
          EPBM0=EPBM
        ENDIF
        RS  = 2*SQRT(EMBM*EPBM)
        IF ( RS.LT.20. ) THEN
          FUNC = 0
          RETURN
        ENDIF
        BMENGY(0) = EBM
        BMENGY(1) = EMBM0
        BMENGY(2) = EPBM0
        BMENGY(3) = EMBM
        BMENGY(4) = EPBM
        WAT = WAT*DPDEBM

      ENDIF
C--
#endif

C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
#if 0
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(19)*1.E6
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
#endif
C--
      QV(0,1)=EMBM+EPBM
      QV(1,1)=0
      QV(2,1)=0
      QV(3,1)=EMBM-EPBM
      IF ( ISRBM.GE.2 ) THEN
        IF( ISR_LLA_ORDER .LE. 0 ) THEN
          CALL ISR_SIMPLE_SPECTRUM(Z(3), Z(4), xALF0,
     >               QED, RS, EMBM, EPBM, QV(0,1))
        ELSE

          CALL ISR_WHIZARD_SPECTRUM(Z(3),Z(4),Z(23),Z(24),Z(25),Z(26),
     >      xALF0, ISR_LLA_ORDER, ROOTS, EMBM, EPBM, QV(0,1),
     >     PISR, WAT, RS)
           IF( RS < 0.0 ) THEN
              FUNC=0
              RETURN
           ENDIF
           DO 10500 KK=0,3
              PVOUT(KK,13)=PISR(KK,1)
              PVOUT(KK,14)=PISR(KK,2)
10500      CONTINUE
        ENDIF
      ENDIF

C--

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
CTTOTSW
CTTOTSW       IDPTCL(1, 3) = 2
CTTOTSW       IDPTCL(2, 3) = 2
CTTOTSW       IDPTCL(3, 3) = 2
C>>>
      IDPTCL(1, 6) = 3
      IDPTCL(2, 6) = 2
      IDPTCL(3, 6) = 2
C
C-- W from t-bar.
      XW1 = XW1MN + Z(20)*(XW1MX-XW1MN)
      DO 200 IMD = IMD1MN, IMD1MX
         IF ( XW1.LE.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   CONTINUE 
#if 0
      IDPTCL(1, 4) = MDVDK(1,IMD,1)
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = MDVDK(3,IMD,1)
      IDPTCL(1, 5) = MDVDK(2,IMD,1)
      IDPTCL(2, 5) = 2
      IDPTCL(3, 5) = MDVDK(3,IMD,1)
#endif
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
310   CONTINUE
#if 0
      IDPTCL(1, 7) = MDVDK(1,IMD,1)
      IDPTCL(2, 7) = 1
      IDPTCL(3, 7) = MDVDK(3,IMD,1)
      IDPTCL(1, 8) = MDVDK(2,IMD,1)
      IDPTCL(2, 8) = 2
      IDPTCL(3, 8) = MDVDK(3,IMD,1)
#endif
C--
      BRW2 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD2 = IMD
      WAT  = WAT*(XW2MX-XW2MN)/BRW2
C--
#ifdef WITH_DBD_STANDARD
      IF( NDCYSL .NE. 0 ) THEN
        XW3 = (( XW2-BRVDK(IMD2-1,1)))/BRW2
        IF( XW3 .GT. 0.5 ) THEN
          LMD1=IMD1
          IMD1=IMD2
          IMD2=LMD1
        ENDIF
        WAT = WAT*2.0D0
      ENDIF
#endif
      IDPTCL(1, 4) = MDVDK(1,IMD1,1)
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = MDVDK(3,IMD1,1)
      IDPTCL(1, 5) = MDVDK(2,IMD1,1)
      IDPTCL(2, 5) = 2
      IDPTCL(3, 5) = MDVDK(3,IMD1,1)
      IDPTCL(1, 7) = MDVDK(1,IMD2,1)
      IDPTCL(2, 7) = 1
      IDPTCL(3, 7) = MDVDK(3,IMD2,1)
      IDPTCL(1, 8) = MDVDK(2,IMD2,1)
      IDPTCL(2, 8) = 2
      IDPTCL(3, 8) = MDVDK(3,IMD2,1)
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( Z(5).LT.HLM ) THEN
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
      Q3MN  = AM3 + AM4 + AM5 + 3*EPS
      Q3MX  = RS - (AM6+AM7+AM8) - 3*EPS
      IF ( Q3MX .LE. Q3MN ) THEN 
        FUNC=0
        RETURN 
      ENDIF
      TH3MN = ATAN((Q3MN-AMT)*(Q3MN+AMT)/AMGT)
      TH3MX = ATAN((Q3MX-AMT)*(Q3MX+AMT)/AMGT)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(6)
      Q32   = AMGT*TAN(TH3) + AMT2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMT2)**2 + AMGT2 )/AMGT
C-- t.
      Q4MN  = AM6 + AM7 + AM8 + 3*EPS
      Q4MX  = RS - Q3 - 3*EPS
      IF ( Q4MX .LE. Q4MN ) THEN
        FUNC=0
        RETURN
      ENDIF
      TH4MN = ATAN((Q4MN-AMT)*(Q4MN+AMT)/AMGT)
      TH4MX = ATAN((Q4MX-AMT)*(Q4MX+AMT)/AMGT)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(7)
      Q42   = AMGT*TAN(TH4) + AMT2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMT2)**2 + AMGT2 )/AMGT
C-- W-.
      Q7MN  = AM4 + AM5 + 2*EPS
      Q7MX  = Q3 - AM3 - 2*EPS
      IF ( Q7MX .LE. Q7MN ) THEN 
        FUNC=0
        RETURN
      ENDIF
      TH7MN = ATAN((Q7MN-AMW)*(Q7MN+AMW)/AMGW)
      TH7MX = ATAN((Q7MX-AMW)*(Q7MX+AMW)/AMGW)
      TH7   = TH7MN + (TH7MX-TH7MN)*Z(8)
      Q72   = AMGW*TAN(TH7) + AMW2
      WAT   = WAT*(TH7MX-TH7MN)*( (Q72-AMW2)**2 + AMGW2 )/AMGW
C-- W+.
      Q8MN  = AM7 + AM8 + 2*EPS
      Q8MX  = Q4 - AM6 - 2*EPS
      IF ( Q8MX .LE. Q8MN ) THEN 
        FUNC=0
        RETURN
      ENDIF 
      TH8MN = ATAN((Q8MN-AMW)*(Q8MN+AMW)/AMGW)
      TH8MX = ATAN((Q8MX-AMW)*(Q8MX+AMW)/AMGW)
      TH8   = TH8MN + (TH8MX-TH8MN)*Z(9)
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
         XPHASE(2,IB) = -1 + 2*Z(2*IB+8)
         XPHASE(3,IB) = x2PI  *Z(2*IB+9)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Calculate differential cross sections.
C--
      CALL SGTTF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
CCDBG       CALL UVZERO(4,PVTEST)
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
C      QV(0,1) = EMBM + EPBM
C      QV(1,1) = 0
C      QV(2,1) = 0
C      QV(3,1) = EMBM - EPBM
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
      COST  = XPHASE(2,1)
      PHIT  = XPHASE(3,1)
      COSBB = XPHASE(2,2)
      PHIBB = XPHASE(3,2)
      COSB  = XPHASE(2,4)
      PHIB  = XPHASE(3,4)
      PTOP  = (RS/2)*BETA(Q42/RS**2,Q32/RS**2)
      CALL PSUMxx(PVOUT(0,3),PVOUT(0, 4),QV(0,3))
      CALL PSUMxx(QV(0,3)   ,PVOUT(0, 5),QV(0,2))
      CALL PSUMxx(PVOUT(0,6),PVOUT(0, 7),QV(0,5))
      CALL PSUMxx(QV(0,5)   ,PVOUT(0, 8),QV(0,4))
      PT1  = SQRT(QV(1,2)**2+QV(2,2)**2+QV(3,2)**2)
      PT2  = SQRT(QV(1,4)**2+QV(2,4)**2+QV(3,4)**2)
C--
      CALL XHFILL( 'h01', DBLE(COST)        ,FUNC )
      CALL XHFILL( 'h02', DBLE(PHIT)*xRD2DG ,FUNC )
      CALL XHFILL( 'h03', DBLE(COSBB)       ,FUNC )
      CALL XHFILL( 'h04', DBLE(PHIBB)*xRD2DG,FUNC )
      CALL XHFILL( 'h05', DBLE(COSB )       ,FUNC )
      CALL XHFILL( 'h06', DBLE(PHIB )*xRD2DG,FUNC )
      CALL XHFILL( 'h07', DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 'h08', DBLE(ROOTS-2*AMT) ,FUNC )
      CALL XHFILL( 'h09', DBLE(Q42)         ,FUNC )
      CALL XHFILL( 'h10', DBLE(Q32)         ,FUNC )
      CALL XHFILL( 'h11', DBLE(PTOP)        ,FUNC )
      CALL XHFILL( 'h12', DBLE(PT1 )        ,FUNC )
      CALL XHFILL( 'h12', DBLE(PT2 )        ,FUNC )
      CALL XHFILL( 'h13', DBLE(IMD1+.1)     ,FUNC )
      CALL XHFILL( 'h14', DBLE(IMD2+.1)     ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
