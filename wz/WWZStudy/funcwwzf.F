C* (Update Record)
C*  95/03/12  K.Fujii		Converted to unix. 
C*				To avoid problems related to kinematic
C*				boundaries, we introduced a tiny offset
C*				EPSL.Now use INSMCN instead of CONINT.

      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "WWZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 8, NHEL = 2 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,5)
      REAL   *8   WAT, WATS
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   -1, +1,   +1,-1,   -1,+1,
     .             -1,+1,   -1, +1,   +1,-1,   +1,-1,
     .             +1,-1,   -1, +1,   +1,-1,   -1,+1,
     .             +1,-1,   -1, +1,   +1,-1,   +1,-1/
      DATA NCALL / 0 /
      EXTERNAL   XHFILL, DHFILL
C--
      DATA EPSL / 1.0 /
C--
#ifdef DEBUG
      INTEGER*4  JDP(2,NP)
#endif
C
C========< Entry Point >================================================
C
#ifdef DEBUG
      print *, '>>>> FUNC >>>>'
#endif
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
         GMH   = GMHTOT
         AMGH  = AMH*GMH
         AMH2  = AMH*AMH
         AMGH2 = AMHZ*AMHZ
C--
         AME   = AMF(1,2,1)
         xALF0 = 1/137.0359895E0
C--
         XW1MN = BRVDK(IMD1MN-1,1)
         XW1MX = BRVDK(IMD1MX  ,1)
         XW2MN = BRVDK(IMD2MN-1,1)
         XW2MX = BRVDK(IMD2MX  ,1)
         XZ0MN = BRVDK(IMD3MN-1,2)
         XZ0MX = BRVDK(IMD3MX  ,2)
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
C      ( 7) : helicity combination for final states
C      ( 8) : m(W-)**2
C      ( 9) : m(W+)**2
C      (10) : m(Z)**2
C      (11) : m(W-W+)**2
C      (12) : cos(theta_Z)
C      (13) : phi_Z
C      (14) : cos(theta_W-)     in W-W+ rest frame
C      (15) : phi_W-            in W-W+ rest frame
C      (16) : cos(theta_fd)     in W-   rest frame
C      (17) : phi_fd            in W-   rest frame
C      (18) : cos(theta_fd_bar) in W+   rest frame
C      (19) : phi_fd_bar        in W+   rest frame
C      (20) : cos(theta_f_bar)  in Z    rest frame
C      (21) : phi_f_bar         in Z    rest frame
C      (22) : final state combination.
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Beam energy.
C--
      RS        = ROOTS
      EBM       = ROOTS/2
      EMBM      = EBM
      EPBM      = EBM
      QED       = 1
C--
C  First decide e- and e+ energy.
C--
      IF ( ISRBM.EQ.3 ) THEN
         call jsfbeamgenr(Z(1), Z(3), Z(2), Z(4), ebm,
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
         DPDEBM    = PEM*PEP
      ELSE
         BMENGY(0) = EBM
         BMENGY(1) = EBM
         BMENGY(2) = EBM
         BMENGY(3) = EBM
         BMENGY(4) = EBM
         DPDEBM    = 1
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
C-- Z.
      XZ = XZ0MN + Z(22)*(XZ0MX-XZ0MN)
      DO 100 IMD = IMD3MN, IMD3MX
         IF ( XZ.LT.BRVDK(IMD,2) )               GO TO 110
100   CONTINUE
C--
110   IDPTCL(1,7) = MDVDK(1,IMD,2)
      IDPTCL(2,7) = MDVDK(2,IMD,2)
      IDPTCL(3,7) = MDVDK(3,IMD,2)
      IDPTCL(1,8) = MDVDK(1,IMD,2)
      IDPTCL(2,8) = MDVDK(2,IMD,2)
      IDPTCL(3,8) = MDVDK(3,IMD,2)
C--
      BRZ  = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD3 = IMD
      WAT  = WAT*(XZ0MX-XZ0MN)/BRZ
C
C-- W-.
      XW1 = XW1MN + ((XZ-BRVDK(IMD-1,2))/BRZ)*(XW1MX-XW1MN)
      DO 200 IMD = IMD1MN, IMD1MX
         IF ( XW1.LT.BRVDK(IMD,1) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,3) = MDVDK(2,IMD,1)
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
      DO 300 IMD = IMD2MN, IMD2MX
         IF ( XW2.LT.BRVDK(IMD,1) )              GO TO 310
300   CONTINUE
C--
310   IDPTCL(1,5) = MDVDK(2,IMD,1)
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
      JHL  = NHEL*Z(7)+1
      JHL  = MIN(JHL,NHEL)
      ICMB = JHL + IOFF
      IF ( IDPTCL(2,7).NE.1 .OR. IDPTCL(3,7).NE.1 ) THEN
         WAT  = WAT*NHEL
      ELSE
         ICMB = 2*((ICMB+1)/2)
         WAT  = WAT*NHEL/2
      ENDIF
      CALL UVCOPY(NP,IHEL(1,ICMB),IHLCMB)
      WATS  = WAT
C--
C  Set invariant masses squared.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AM7   = AMF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AM8   = AMF(IDPTCL(1,8),IDPTCL(2,8),IDPTCL(3,8))
C-- W-.
      Q4MN  = AM3 + AM4 + EPSL
      Q4MX  = MAX( RS - (AM5+AM6+AM7+AM8) - EPSL, Q4MN )
      TH4MN = ATAN((Q4MN-AMW)*(Q4MN+AMW)/AMGW)
      TH4MX = ATAN((Q4MX-AMW)*(Q4MX+AMW)/AMGW)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(8)
      Q42   = AMGW*TAN(TH4) + AMW2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMW2)**2 + AMGW2 )/AMGW
C-- W+.
      Q5MN  = AM5 + AM6 + EPSL
      Q5MX  = MAX( RS - (Q4+AM7+AM8) - EPSL, Q5MN )
      TH5MN = ATAN((Q5MN-AMW)*(Q5MN+AMW)/AMGW)
      TH5MX = ATAN((Q5MX-AMW)*(Q5MX+AMW)/AMGW)
      TH5   = TH5MN + (TH5MX-TH5MN)*Z(9)
      Q52   = AMGW*TAN(TH5) + AMW2
      Q5    = SQRT(Q52)
      WAT   = WAT*(TH5MX-TH5MN)*( (Q52-AMW2)**2 + AMGW2 )/AMGW
C-- Z.
      Q2MN  = AM7 + AM8 + EPSL
      Q2MX  = MAX( RS - (Q4+Q5) - EPSL, Q2MN )
      TH2MN = ATAN((Q2MN-AMZ)*(Q2MN+AMZ)/AMGZ)
      TH2MX = ATAN((Q2MX-AMZ)*(Q2MX+AMZ)/AMGZ)
      TH2   = TH2MN + (TH2MX-TH2MN)*Z(10)
      Q22   = AMGZ*TAN(TH2) + AMZ2
      Q2    = SQRT(Q22)
      WAT   = WAT*(TH2MX-TH2MN)*( (Q22-AMZ2)**2 + AMGZ2 )/AMGZ
C-- W-W+.
#ifdef HIGGSONLY
      Q32MN = Q4+Q5
      Q32MX = MAX( RS-Q2, Q32MN )
      TH32MN = ATAN((Q32MN-AMH)*(Q32MN+AMH)/AMGH)
      TH32MX = ATAN((Q32MX-AMH)*(Q32MX+AMH)/AMGH)
      TH32   = TH32MN + (TH32MX-TH32MN)*Z(11)
      Q32    = AMGH*TAN(TH32) + AMH2
      WAT   = WAT*(TH32MX-TH32MN)*( (Q32-AMH2)**2 + AMGH2 )/AMGH
#else
      Q32MN = (Q4+Q5)**2
      Q32MX = MAX( (RS-Q2)**2, Q32MN )
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(11)
      WAT   = WAT*(Q32MX-Q32MN)
#endif
C--
#ifdef ZEROWIDTH
      Q32MN = (AMW+AMW)**2
      Q32MX = (RS-AMZ)**2
      Q22   = AMZ2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(11)
      Q42   = AMW2
      Q52   = AMW2
C--
      WAT   = WATS*(Q32MX-Q32MN)
     .            *(xPI*AMW*GMW)**2
     .            *(xPI*AMZ*GMZ)
#endif
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q32
      XPHASE(1,3) = Q42
      XPHASE(1,4) = Q52
      XPHASE(1,5) = Q22
C--
C  Angular variables.
C--
      DO 400 IB = 1, 5
         XPHASE(2,IB) = -1 + 2*Z(2*IB+10)
         XPHASE(3,IB) = x2PI  *Z(2*IB+11)
         WAT = WAT*2*x2PI
400   CONTINUE
#ifdef DEBUG
      DO 665 K = 1, NP
         JDP(1,K) = 100*IDPTCL(1,K)+10*IDPTCL(2,K)+IDPTCL(3,K)
         JDP(2,K) = IHLCMB(K)
665   CONTINUE
      PRINT *, ' -------    '
      PRINT *, ' IDP, IHL = ', ((JDP(K,KK),K=1,2),KK=1,NP)
      DO 666 KK = 1, NP-3
         PRINT *, ' XPHASE(*,', KK, ') = ', (XPHASE(K,KK),K=1,3)
666   CONTINUE
#endif
C--
C  Calculate differential cross sections.
C--
      CALL SGWWZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG     = SG*QED
      FUNC   = SG*WAT
      DSGMDX = SG
#ifdef DEBUG
      PRINT *, ' SG = ', SG,' FUNC = ', FUNC
#endif
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
      IF ( ABS(QV(3,1)).GT.1.E-3 ) THEN
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
      CALL PSUMxx(PV(0,3),PV(0,4),QV(0,4))
      CALL PSUMxx(PV(0,5),PV(0,6),QV(0,5))
      CALL PSUMxx(PV(0,7),PV(0,8),QV(0,2))
      AMWW2 = XPHASE(1,2)
      AMWW  = SQRT(AMWW2)
      EBM   = ROOTS/2
      EWM   = QV(0,4)
      EWP   = QV(0,5)
      EZ    = QV(0,2)
      RDDG  = 180/xPI
      COSZ  = XPHASE(2,1)
      PHIZ  = XPHASE(3,1)
      COSWM = XPHASE(2,2)
      PHIWM = XPHASE(3,2)
      CALL XHFILL( 'h01', DBLE(COSZ)        ,FUNC )
      CALL XHFILL( 'h02', DBLE(PHIZ*RDDG)   ,FUNC )
      CALL XHFILL( 'h03', DBLE(AMWW/ROOTS)  ,FUNC )
      CALL XHFILL( 'h04', DBLE(COSWM)       ,FUNC )
      CALL XHFILL( 'h05', DBLE(PHIWM*RDDG)  ,FUNC )
      CALL XHFILL( 'h07', DBLE(ICMB)        ,FUNC )
      CALL DHFILL( 'hd09', DBLE(EZ/EBM),  DBLE(COSZ)   , FUNC )
      CALL DHFILL( 'hd10', DBLE(EWM/EBM), DBLE(EWP/EBM), FUNC )
C--
C  That's it.
C--
      RETURN
      END
