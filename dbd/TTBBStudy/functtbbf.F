C* 2009/07/30  T.Tanabe   FUNC for e+ e- --> t t_bar g --> t t_bar b b_bar full amplitudes.
C*                        Based on TTHStudy by K.Fujii
C*

      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
#include "TTBBCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
#include "USRPRM.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 10, NHEL = 2 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,10)
C--
      REAL   *8   WAT, WATS
      REAL   *8   Q3MN, Q3MX, TH3MN, TH3MX, TH3,
     .            Q4MN, Q4MX, TH4MN, TH4MX, TH4,
     .            Q5MN, Q5MX, TH5MN, TH5MX, TH5,
     .            Q7MN, Q7MX, TH7MN, TH7MX, TH7,
     .            Q8MN, Q8MX, TH8MN, TH8MX, TH8
      REAL   *4   M2TG, M2TBG
#ifdef WITH_DBD_STANDARD
      REAL*8 PISR(0:3,2)
#endif
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   +1, +1,-1,   -1, -1,+1,   -1,+1,
     .             -1,+1,   +1, +1,-1,   -1, -1,+1,   +1,-1,
     .             +1,-1,   +1, +1,-1,   -1, -1,+1,   -1,+1,
     .             +1,-1,   +1, +1,-1,   -1, -1,+1,   +1,-1/
#ifdef DEBUG      
      INTEGER*4  JDP(2,NP)
#endif
      DATA NCALL / 0 /
      EXTERNAL   XHFILL, DHFILL
C--
c      DATA EPSL  / 0.0 /
      DATA EPSL  / 5.D0 /
C
C========< Entry Point >================================================
C
#ifdef DEBUG
       print *, '>>>> Now enters FUNC.'
#endif
C--
C  Set some variables.
C     ENGM  = variable mass range in units of width.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         AMT   = AMF(3,1,2)
         GMT   = GMF(3,1,2)
         AMGT  = AMT*GMT
         AMT2  = AMT*AMT
         AMGT2 = AMGT*AMGT
C--
         GMW   = GMWTOT
         AMGW  = AMW*GMW
         AMW2  = AMW*AMW
         AMGW2 = AMGW*AMGW
C--
         GMH   = GMHTOT
         AMGH  = AMH*GMH
         AMH2  = AMH*AMH
         AMGH2 = AMGH*AMGH
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
C     Z( 1) : e- beam
C      ( 2) : e+ beam
C      ( 3) : ISR for e-
C      ( 4) : ISR for e+ ( or e-/e+ selection when ISR_LLA_Order=0 )
C      ( 5) : e- helicity
C      ( 6) : helicity combination for final states.
C      ( 7) : m(t-bar)**2
C      ( 8) : m(t)**2
C      ( 9) : m(g)**2
C      (10) : m(t-t_bar)**2
C      (11) : m(W-)**2
C      (12) : m(W+)**2
C      (13) : cos(theta_g)
C      (14) : phi_H
C      (15) : cos(theta_t)     in t-t_bar rest frame
C      (16) : phi_t            in t-t_bar rest frame
C      (17) : cos(theta_q_bar) in t_bar rest frame
C      (18) : phi_q_bar        in t_bar rest frame
C      (19) : cos(theta_f)     in W- rest frame
C      (20) : phi_f            in W- rest frame
C      (21) : cos(theta_q)     in t rest frame
C      (22) : phi_q            in t rest frame
C      (23) : cos(theta_f_bar) in W+ rest frame
C      (24) : phi_f_bar        in W+ rest frame
C      (25) : cos(theta_f_bar) in Z rest frame
C      (26) : phi_f_bar        in Z rest frame
C      (27) : final state combination.
C      (28) : e- beam gaussian spread
C      (29) : e+ beam gaussian spread
C      (30) : random variable for Pt of ISR photon from e- beam
C      (31) : random variable for Pt of ISR photon from e+ beam
C      (32) : random variable for azimuthal angle of ISR photon from e- beam
C      (33) : random variable for azimuthal angle of ISR photon from e+ beam
C--
      WAT = 1
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
#ifndef WITH_DBD_STANDARD
      IF ( ISRBM.EQ.3 ) THEN
#if 0
         CALL PROBEB(REAL(Z(1)),REAL(Z(27)),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(REAL(Z(2)),REAL(Z(28)),EBM,EPBM,EPBM0,PEP)
#else
         call jsfbeamgenr(z(1), z(2), z(27), z(28), ebm,
     >        embm, epbm, embm0, epbm0, dpdebm )
#endif
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
#if 0
         DPDEBM    = PEM*PEP
#else
         wat = wat*dpdebm
#endif
      ENDIF
#else
C -- when WITH_DBD_STANDARD is set
      IF ( ISRBM.GE.3 ) THEN
        IF ( ISRBM.EQ.3 ) THEN
          CALL JSFBEAMGENR(Z(1), Z(2), Z(28), Z(29), EBM,
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
C  Then decide reduced sqrt(s) after ISR.
C--
      QED = 1
C--
#if 0
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(26)*1.E6
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
          CALL ISR_WHIZARD_SPECTRUM(Z(3),Z(4),Z(30),Z(31),Z(32),Z(33),
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
C     t --> b + W 
C     g --> b + bbar only.
C--
      IDPTCL(1, 3) = 3
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 2
      IDPTCL(1, 6) = 3
      IDPTCL(2, 6) = 2
      IDPTCL(3, 6) = 2
C-- 
      IDPTCL(1, 9) = 3
      IDPTCL(2, 9) = 2
      IDPTCL(3, 9) = 2
      IDPTCL(1,10) = 3
      IDPTCL(2,10) = 2
      IDPTCL(3,10) = 2
C
C-- W from t-bar.
      XW1 = XW1MN + Z(27)*(XW1MX-XW1MN)
      DO 200 IMD = IMD1MN, IMD1MX
         IF ( XW1.LE.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   CONTINUE
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
C--
      BRW2 = BRVDK(IMD,1) - BRVDK(IMD-1,1)
      IMD2 = IMD
      WAT  = WAT*(XW2MX-XW2MN)/BRW2
C
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
C  Select electron helicity.
C--
      HLM   = (1-POLE)/2
      IF ( Z(5).LT.HLM ) THEN
         IOFF = 0
         WAT = WAT*(1+POLP)/2
      ELSE
         IOFF = NHEL
         WAT = WAT*(1-POLP)/2
      ENDIF
C--
C  Select final state helicities.
C--
      JHL  = NHEL*Z(6)+1
      JHL  = MIN(JHL,NHEL)
      ICMB = JHL + IOFF
      WAT  = WAT*NHEL
      CALL UVCOPY(NP,IHEL(1,ICMB),IHLCMB)
      WATS  = WAT
#ifdef DEBUG
       print *, 'Now calculate kinematics'
#endif
C--
C  Set .
C--
      AM3   = AMF(IDPTCL(1, 3),IDPTCL(2, 3),IDPTCL(3, 3))
      AM4   = AMF(IDPTCL(1, 4),IDPTCL(2, 4),IDPTCL(3, 4))
      AM5   = AMF(IDPTCL(1, 5),IDPTCL(2, 5),IDPTCL(3, 5))
      AM6   = AMF(IDPTCL(1, 6),IDPTCL(2, 6),IDPTCL(3, 6))
      AM7   = AMF(IDPTCL(1, 7),IDPTCL(2, 7),IDPTCL(3, 7))
      AM8   = AMF(IDPTCL(1, 8),IDPTCL(2, 8),IDPTCL(3, 8))
      AM9   = AMF(IDPTCL(1, 9),IDPTCL(2, 9),IDPTCL(3, 9))
      AM10  = AMF(IDPTCL(1,10),IDPTCL(2,10),IDPTCL(3,10))
C--
C  Check if RS exceeds the sum of the masses.
C--
      AMSUM = AM3 + AM4 + AM5 + AM6 + AM7 + AM8 + AM9 + AM10
      IF ( RS.LE.AMSUM+EPSL ) THEN
         FUNC = 0
         RETURN
      ENDIF
C--
C  Set invariant masses squared.
C--
#ifdef DEBUG
       print *, '--tbar'
#endif
C-- t_bar.
      Q3MN  = AM3 + AM4 + AM5 + EPSL
      Q3MX  = MAX( RS - (AM6+AM7+AM8+AM9+AM10) - EPSL, Q3MN )
      TH3MN = ATAN((Q3MN-AMT)*(Q3MN+AMT)/AMGT)
      TH3MX = ATAN((Q3MX-AMT)*(Q3MX+AMT)/AMGT)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(7)
      Q32   = AMGT*TAN(TH3) + AMT2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMT2)**2 + AMGT2 )/AMGT
#ifdef DEBUG
       print *, '--t'
#endif
C-- t.
      Q4MN  = AM6 + AM7 + AM8 + EPSL
      Q4MX  = MAX( RS - (Q3+AM9+AM10) - EPSL, Q4MN )
      TH4MN = ATAN((Q4MN-AMT)*(Q4MN+AMT)/AMGT)
      TH4MX = ATAN((Q4MX-AMT)*(Q4MX+AMT)/AMGT)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(8)
      Q42   = AMGT*TAN(TH4) + AMT2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMT2)**2 + AMGT2 )/AMGT
#ifdef DEBUG
       print *, '--g'
#endif
C-- g.
      Q5MN  = AM9 + AM10 + EPSL
      Q5MX  = MAX( RS - (Q3+Q4) - EPSL, Q5MN )
      Q52MN = Q5MN**2
      Q52MX = Q5MX**2
      Q52   = Q52MN + (Q52MX-Q52MN)*Z(9)
      Q5    = SQRT(Q52)
      WAT   = WAT*(Q52MX-Q52MN)
#ifdef DEBUG
       print *, '--ttbar'
#endif
C-- t-t_bar.
      Q62MN = (Q3+Q4)**2
      Q62MX = (RS-Q5)**2
      Q62   = Q62MN + (Q62MX-Q62MN)*Z(10)
      WAT   = WAT*(Q62MX-Q62MN)
#ifdef DEBUG
       print *, '--W-'
#endif
C-- W-.
      Q7MN  = AM4 + AM5 + EPSL
      Q7MX  = MAX( Q3 - AM3 - EPSL, Q7MN )
      TH7MN = ATAN((Q7MN-AMW)*(Q7MN+AMW)/AMGW)
      TH7MX = ATAN((Q7MX-AMW)*(Q7MX+AMW)/AMGW)
      TH7   = TH7MN + (TH7MX-TH7MN)*Z(11)
      Q72   = AMGW*TAN(TH7) + AMW2
      WAT   = WAT*(TH7MX-TH7MN)*( (Q72-AMW2)**2 + AMGW2 )/AMGW
#ifdef DEBUG
       print *, '--W+'
#endif
C-- W+.
      Q8MN  = AM7 + AM8 + EPSL
      Q8MX  = MAX( Q4 - AM6 - EPSL, Q8MN )
      TH8MN = ATAN((Q8MN-AMW)*(Q8MN+AMW)/AMGW)
      TH8MX = ATAN((Q8MX-AMW)*(Q8MX+AMW)/AMGW)
      TH8   = TH8MN + (TH8MX-TH8MN)*Z(12)
      Q82   = AMGW*TAN(TH8) + AMW2
      WAT   = WAT*(TH8MX-TH8MN)*( (Q82-AMW2)**2 + AMGW2 )/AMGW
#ifdef DEBUG
       print *, 'WAT = ', WAT
#endif
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q62
      XPHASE(1,3) = Q32
      XPHASE(1,4) = Q72
      XPHASE(1,5) = Q42
      XPHASE(1,6) = Q82
      XPHASE(1,7) = Q52
C--
C  Angular variables.
C--
      DO 400 IB = 1, 7
         XPHASE(2,IB) = -1 + 2*Z(2*IB+10)
         XPHASE(3,IB) = x2PI  *Z(2*IB+11)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Calculate differential cross sections.
C--
#ifdef DEBUG
       print *, 'Going to call SGTTBBF'
#endif
      CALL SGTTBBF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG   = SG*QED
      FUNC = SG*WAT
#ifdef DEBUG 
      DO 665 K = 1, NP
         JDP(1,K) = 100*IDPTCL(1,K)+10*IDPTCL(2,K)+IDPTCL(3,K)
         JDP(2,K) = IHLCMB(K)
665   CONTINUE
      PRINT *, ' -------    '
      PRINT *, ' IDP, IHL = ', ((JDP(K,KK),K=1,2),KK=1,NP)
      DO 666 KK = 1, 7
         PRINT *, ' XPHASE(*,', KK, ') = ', (XPHASE(K,KK),K=1,3)
666   CONTINUE
      PRINT *, ' SG = ', SG,' FUNC = ', FUNC, ' QED = ', QED
#endif
      IF ( SG.LE.0. ) RETURN
C--
C  Boost everything to laboratroy frame.
C--
      IF ( ABS(QV(3,1)).GT.1.E-3 ) THEN
#ifdef DEBUG
      PRINT *, ' QV = ', (QV(K,1),K=0,3)
#endif
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
#ifdef DEBUG
       print *, 'Hists and plots.'
#endif
      CALL PSUMxx(PV(0, 3),PV(0, 4),QV(0,3))
      CALL PSUMxx(QV(0, 3),PV(0, 5),QV(0,3))
      CALL PSUMxx(PV(0, 6),PV(0, 7),QV(0,4))
      CALL PSUMxx(QV(0, 4),PV(0, 8),QV(0,4))
      CALL PSUMxx(PV(0, 9),PV(0,10),QV(0,5))
C---  tbar + g
      CALL PSUMxx(QV(0,3),QV(0,5),QV(0,6))
      CALL RDOTxx(QV(0,6),QV(0,6),M2TBG);
C---  t + g
      CALL PSUMxx(QV(0,4),QV(0,5),QV(0,7))
      CALL RDOTxx(QV(0,7),QV(0,7),M2TG);
C--
      AMSTT2 = XPHASE(1,2)
      AMTT  = SQRT(AMSTT2)
      ETB   = QV(0,3)
      ET    = QV(0,4)
      EZ    = QV(0,5)
      RDDG  = 180/xPI
      APE   = SQRT( PV(1,1)**2 + PV(2,1)**2 + PV(3,1)**2 )
      APT   = SQRT( QV(1,4)**2 + QV(2,4)**2 + QV(3,4)**2 )
      CST   = (QV(1,4)*PV(1,1)+QV(2,4)*PV(2,1)+QV(3,4)*PV(3,1))/APT/APE
      COSZ  = XPHASE(2,1)
      PHIZ  = XPHASE(3,1)
      COST  = XPHASE(2,2)
      PHIT  = XPHASE(3,2)
      QMWM  = SQRT(Q72)
      QMWP  = SQRT(Q82)
      CALL XHFILL( 'h01', DBLE(COSZ)        ,FUNC )
      CALL XHFILL( 'h02', DBLE(PHIZ)*RDDG   ,FUNC )
#if 0
      CALL XHFILL( 'h03', DBLE(AMTT/ROOTS)  ,FUNC )
#else
      CALL XHFILL( 'h03', DBLE(AMTT)        ,FUNC )
#endif
      CALL XHFILL( 'h04', DBLE(COST)        ,FUNC )
      CALL XHFILL( 'h05', DBLE(PHIT)*RDDG   ,FUNC )
      CALL XHFILL( 'h06', DBLE(CST)         ,FUNC )
      CALL XHFILL( 'h07', DBLE(ICMB)        ,FUNC )
      CALL XHFILL( 'h08', DBLE(Q3)          ,FUNC )
      CALL XHFILL( 'h09', DBLE(Q4)          ,FUNC )
      CALL XHFILL( 'h10', DBLE(Q5)          ,FUNC )
      CALL XHFILL( 'h11', DBLE(QMWM)        ,FUNC )
      CALL XHFILL( 'h12', DBLE(QMWP)        ,FUNC )
      CALL XHFILL( 'h13', DBLE(RS/ROOTS)    ,FUNC )
      CALL DHFILL( 'hd21', DBLE(EZ/EBM), DBLE(COSZ)    , FUNC )
      CALL DHFILL( 'hd22', DBLE(ET/EBM), DBLE(ETB/EBM) , FUNC )
      CALL XHFILL( 'h15', DBLE(IMD1+.1)     ,FUNC )
      CALL XHFILL( 'h16', DBLE(IMD2+.1)     ,FUNC )
      CALL DHFILL( 'dalitz', DBLE(M2TG), DBLE(M2TBG) , FUNC )
#ifdef DEBUG
       print *, 'That''s it.'
#endif
C--
C  That's it.
C--
      RETURN
      END
