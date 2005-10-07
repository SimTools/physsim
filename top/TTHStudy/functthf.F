C* 95/03/28  K.Fujii    FUNC for e+ e- --> t t_bar H full amplitudes.
C* 95/06/07  K.Fujii	Use ISRBM switch.
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
      PARAMETER  ( NP = 10, NHEL = 2 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,5)
C--
      REAL   *8   WAT, WATS
      REAL   *8   Q3MN, Q3MX, TH3MN, TH3MX, TH3,
     .            Q4MN, Q4MX, TH4MN, TH4MX, TH4,
     .            Q5MN, Q5MX, TH5MN, TH5MX, TH5,
     .            Q7MN, Q7MX, TH7MN, TH7MX, TH7,
     .            Q8MN, Q8MX, TH8MN, TH8MX, TH8
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   +1, +1,-1,   -1, -1,+1,   -1,-1,
     .             -1,+1,   +1, +1,-1,   -1, -1,+1,   +1,+1,
     .             +1,-1,   +1, +1,-1,   -1, -1,+1,   -1,-1,
     .             +1,-1,   +1, +1,-1,   -1, -1,+1,   +1,+1/
#ifdef DEBUG      
      INTEGER*4  JDP(2,NP)
#endif
      DATA NCALL / 0 /
      EXTERNAL   XHFILL, DHFILL
C--
      DATA EPSL  / 0.0 /
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
C      ( 6) : m(t-bar)**2
C      ( 7) : m(t)**2
C      ( 8) : m(J)**2
C      ( 9) : m(t-t_bar)**2
C      (10) : m(W-)**2
C      (11) : m(W+)**2
C      (12) : cos(theta_H)
C      (13) : phi_H
C      (14) : cos(theta_t)     in t-t_bar rest frame
C      (15) : phi_t            in t-t_bar rest frame
C      (16) : cos(theta_q_bar) in t_bar rest frame
C      (17) : phi_q_bar        in t_bar rest frame
C      (18) : cos(theta_f)     in W- rest frame
C      (19) : phi_f            in W- rest frame
C      (20) : cos(theta_q)     in t rest frame
C      (21) : phi_q            in t rest frame
C      (22) : cos(theta_f_bar) in W+ rest frame
C      (23) : phi_f_bar        in W+ rest frame
C      (24) : cos(theta_f_bar) in Z rest frame
C      (25) : phi_f_bar        in Z rest frame
C      (26) : final state combination.
C      (27) : e- beam gaussian spread
C      (28) : e+ beam gaussian spread
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
      IF ( ISRBM.EQ.3 ) THEN
         CALL PROBEB(REAL(Z(1)),REAL(Z(27)),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(REAL(Z(2)),REAL(Z(28)),EBM,EPBM,EPBM0,PEP)
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
      QED = 1
C--
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
C     H --> b + bbar only.
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
      XW1 = Z(26)
      DO 200 IMD = 1, 12
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 210
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
C
C-- W from t.
      XW2 = (XW1-BRVDK(IMD-1,1))/BRW1
      DO 300 IMD = 1, 12
         IF ( XW2.LT.BRVDK(IMD,1) )              GO TO 310
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
      WAT  = WAT/BRW1/BRW2
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
      IF ( IDPTCL(2,9).NE.1 .OR. IDPTCL(3,9).NE.1 ) THEN
         WAT  = WAT*NHEL
      ELSE
         ICMB = 2*((ICMB+1)/2)
         WAT  = WAT*NHEL/2
      ENDIF
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
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(6)
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
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(7)
      Q42   = AMGT*TAN(TH4) + AMT2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMT2)**2 + AMGT2 )/AMGT
#ifdef DEBUG
       print *, '--Z'
#endif
C-- H.
      Q5MN  = AM9 + AM10 + EPSL
      Q5MX  = MAX( RS - (Q3+Q4) - EPSL, Q5MN )
      TH5MN = ATAN((Q5MN-AMH)*(Q5MN+AMH)/AMGH)
      TH5MX = ATAN((Q5MX-AMH)*(Q5MX+AMH)/AMGH)
      TH5   = TH5MN + (TH5MX-TH5MN)*Z(8)
      Q52   = AMGH*TAN(TH5) + AMH2
      Q5    = SQRT(Q52)
      WAT   = WAT*(TH5MX-TH5MN)*( (Q52-AMH2)**2 + AMGH2 )/AMGH
#ifdef DEBUG
       print *, '--ttbar'
#endif
C-- t-t_bar.
      Q62MN = (Q3+Q4)**2
      Q62MX = (RS-Q5)**2
      Q62   = Q62MN + (Q62MX-Q62MN)*Z(9)
      WAT   = WAT*(Q62MX-Q62MN)
#ifdef DEBUG
       print *, '--W-'
#endif
C-- W-.
      Q7MN  = AM4 + AM5 + EPSL
      Q7MX  = MAX( Q3 - AM3 - EPSL, Q7MN )
      TH7MN = ATAN((Q7MN-AMW)*(Q7MN+AMW)/AMGW)
      TH7MX = ATAN((Q7MX-AMW)*(Q7MX+AMW)/AMGW)
      TH7   = TH7MN + (TH7MX-TH7MN)*Z(10)
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
      TH8   = TH8MN + (TH8MX-TH8MN)*Z(11)
      Q82   = AMGW*TAN(TH8) + AMW2
      WAT   = WAT*(TH8MX-TH8MN)*( (Q82-AMW2)**2 + AMGW2 )/AMGW
#ifdef DEBUG
       print *, 'WAT = ', WAT
#endif
C--
#ifdef ZEROWIDTH
      Q62MN = (AMT+AMT)**2
      Q62MX = (RS-AMH)**2
      Q32   = AMT2
      Q42   = AMT2
      Q52   = AMH2
      Q62   = Q62MN + (Q62MX-Q62MN)*Z(9)
      Q72   = AMW2
      Q82   = AMW2
C--
      WAT   = WATS*(Q62MX-Q62MN)
     .            *(xPI*AMT*GMT)**2
     .            *(xPI*AMW*GMW)**2
     .            *(xPI*AMH*GMH)
#endif ZEROWIDTH
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
       print *, 'Going to call SGTTHF'
#endif
      CALL SGTTHF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
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
      CALL XHFILL( 'h03', DBLE(AMTT/ROOTS)  ,FUNC )
      CALL XHFILL( 'h04', DBLE(COST)        ,FUNC )
      CALL XHFILL( 'h05', DBLE(PHIT)*RDDG   ,FUNC )
      CALL XHFILL( 'h06', DBLE(CST)         ,FUNC )
      CALL XHFILL( 'h07', DBLE(ICMB)        ,FUNC )
      CALL XHFILL( 'h08', DBLE(Q3)          ,FUNC )
      CALL XHFILL( 'h09', DBLE(Q4)          ,FUNC )
      CALL XHFILL( 'h10', DBLE(Q5)          ,FUNC )
      CALL XHFILL( 'h11', DBLE(QMWM)        ,FUNC )
      CALL XHFILL( 'h12', DBLE(QMWP)        ,FUNC )
      CALL DHFILL( 'hd21', DBLE(EZ/EBM), DBLE(COSZ)    , FUNC )
      CALL DHFILL( 'hd22', DBLE(ET/EBM), DBLE(ETB/EBM) , FUNC )
#ifdef DEBUG
       print *, 'That''s it.'
#endif
C--
C  That's it.
C--
      RETURN
      END


