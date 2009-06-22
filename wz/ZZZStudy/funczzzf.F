C* (Update Record)
C*  09/06/18  K.Fujii		Original version derived from WWZ

      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "ZZZCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 8, NHEL = 8 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,5)
      REAL   *8   WAT, WATS
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   -1, +1,   -1,+1,   -1,+1,
     .             -1,+1,   -1, +1,   -1,+1,   +1,-1,
     .             -1,+1,   -1, +1,   +1,-1,   -1,+1,
     .             -1,+1,   -1, +1,   +1,-1,   +1,-1,
     .             -1,+1,   +1, -1,   -1,+1,   -1,+1,
     .             -1,+1,   +1, -1,   -1,+1,   +1,-1,
     .             -1,+1,   +1, -1,   +1,-1,   -1,+1,
     .             -1,+1,   +1, -1,   +1,-1,   +1,-1,
     .             +1,-1,   -1, +1,   -1,+1,   -1,+1,
     .             +1,-1,   -1, +1,   -1,+1,   +1,-1,
     .             +1,-1,   -1, +1,   +1,-1,   -1,+1,
     .             +1,-1,   -1, +1,   +1,-1,   +1,-1,
     .             +1,-1,   +1, -1,   -1,+1,   -1,+1,
     .             +1,-1,   +1, -1,   -1,+1,   +1,-1,
     .             +1,-1,   +1, -1,   +1,-1,   -1,+1,
     .             +1,-1,   +1, -1,   +1,-1,   +1,-1/
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
         XZ1MN = BRVDK(IMD1MN-1,2)
         XZ1MX = BRVDK(IMD1MX  ,2)
         XZ2MN = BRVDK(IMD2MN-1,2)
         XZ2MX = BRVDK(IMD2MX  ,2)
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
C      ( 8) : m(Z1)**2
C      ( 9) : m(Z2)**2
C      (10) : m(Z0)**2
C      (11) : m(Z1Z2)**2
C      (12) : cos(theta_Z0)
C      (13) : phi_Z0
C      (14) : cos(theta_Z1)     in Z1Z2 rest frame
C      (15) : phi_Z1            in Z1Z2 rest frame
C      (16) : cos(theta_fd)     in Z1   rest frame
C      (17) : phi_fd            in Z1   rest frame
C      (18) : cos(theta_fd_bar) in Z2   rest frame
C      (19) : phi_fd_bar        in Z2   rest frame
C      (20) : cos(theta_f_bar)  in Z0   rest frame
C      (21) : phi_f_bar         in Z0   rest frame
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
         wat = wat*dpdebm
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
         IF ( XZ.LE.BRVDK(IMD,2) )               GO TO 110
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
C-- Z1.
      XZ1 = XZ1MN + ((XZ-BRVDK(IMD-1,2))/BRZ)*(XZ1MX-XZ1MN)
      DO 200 IMD = IMD1MN, IMD1MX
         IF ( XZ1.LE.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,3) = MDVDK(1,IMD,2)
      IDPTCL(2,3) = MDVDK(2,IMD,2)
      IDPTCL(3,3) = MDVDK(3,IMD,2)
      IDPTCL(1,4) = MDVDK(1,IMD,2)
      IDPTCL(2,4) = MDVDK(2,IMD,2)
      IDPTCL(3,4) = MDVDK(3,IMD,2)
C--
      BRZ1 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD1 = IMD
      WAT  = WAT*(XZ1MX-XZ1MN)/BRZ1
C
C-- Z2.
      XZ2 = XZ2MN + ((XZ1-BRVDK(IMD-1,2))/BRZ1)*(XZ2MX-XZ2MN)
      DO 300 IMD = IMD2MN, IMD2MX
         IF ( XZ2.LE.BRVDK(IMD,2) )              GO TO 310
300   CONTINUE
C--
310   IDPTCL(1,5) = MDVDK(1,IMD,2)
      IDPTCL(2,5) = MDVDK(2,IMD,2)
      IDPTCL(3,5) = MDVDK(3,IMD,2)
      IDPTCL(1,6) = MDVDK(1,IMD,2)
      IDPTCL(2,6) = MDVDK(2,IMD,2)
      IDPTCL(3,6) = MDVDK(3,IMD,2)
C--
      BRZ2 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD2 = IMD
      WAT  = WAT*(XZ2MX-XZ2MN)/BRZ2
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
      ICMB = JHL
#if 0
C--
C  Avoid right-handed helicity for nu to save CPU time
C  since coupling is zero.
C--
      IF ( IDPTCL(2,3).NE.1 .OR. IDPTCL(3,3).NE.1 ) THEN
         IF ( IDPTCL(2,5).NE.1 .OR. IDPTCL(3,5).NE.1 ) THEN
            IF ( IDPTCL(2,7).NE.1 .OR. IDPTCL(3,7).NE.1 ) THEN
               WAT  = WAT*NHEL
            ELSE
               ICMB = 2*((ICMB+1)/2)
               WAT  = WAT*NHEL/2
            ENDIF
         ELSE
            IF ( IDPTCL(2,7).NE.1 .OR. IDPTCL(3,7).NE.1 ) THEN
               ICMB = 4*((ICMB-1)/4) + MOD(ICMB,2) + 3
               WAT  = WAT*NHEL/2
            ELSE
               ICMB = 4*((ICMB+3)/4)
               WAT  = WAT*NHEL/4
            ENDIF
         ENDIF
      ELSE
         IF ( IDPTCL(2,5).NE.1 .OR. IDPTCL(3,5).NE.1 ) THEN
            IF ( IDPTCL(2,7).NE.1 .OR. IDPTCL(3,7).NE.1 ) THEN
               ICMB = ICMB - 4*((ICMB-1)/4)
               WAT  = WAT*NHEL/2
            ELSE
               ICMB = 2*((ICMB - 4*((ICMB-1)/4) + 1)/2)
               WAT  = WAT*NHEL/4
            ENDIF
         ELSE
            IF ( IDPTCL(2,7).NE.1 .OR. IDPTCL(3,7).NE.1 ) THEN
               ICMB = MOD(ICMB,2) + 3
               WAT  = WAT*NHEL/4
            ELSE
               ICMB = 4
               WAT  = WAT*NHEL/8
            ENDIF
         ENDIF
      ENDIF
#else
      WAT  = WAT*NHEL
#endif
      ICMB = ICMB + IOFF
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
C-- Z1.
      Q4MN  = AM3 + AM4 + EPSL
      Q4MX  = MAX( RS - (AM5+AM6+AM7+AM8) - EPSL, Q4MN )
      TH4MN = ATAN((Q4MN-AMZ)*(Q4MN+AMZ)/AMGZ)
      TH4MX = ATAN((Q4MX-AMZ)*(Q4MX+AMZ)/AMGZ)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(8)
      Q42   = AMGZ*TAN(TH4) + AMZ2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Z2.
      Q5MN  = AM5 + AM6 + EPSL
      Q5MX  = MAX( RS - (Q4+AM7+AM8) - EPSL, Q5MN )
      TH5MN = ATAN((Q5MN-AMZ)*(Q5MN+AMZ)/AMGZ)
      TH5MX = ATAN((Q5MX-AMZ)*(Q5MX+AMZ)/AMGZ)
      TH5   = TH5MN + (TH5MX-TH5MN)*Z(9)
      Q52   = AMGZ*TAN(TH5) + AMZ2
      Q5    = SQRT(Q52)
      WAT   = WAT*(TH5MX-TH5MN)*( (Q52-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Z.
      Q2MN  = AM7 + AM8 + EPSL
      Q2MX  = MAX( RS - (Q4+Q5) - EPSL, Q2MN )
      TH2MN = ATAN((Q2MN-AMZ)*(Q2MN+AMZ)/AMGZ)
      TH2MX = ATAN((Q2MX-AMZ)*(Q2MX+AMZ)/AMGZ)
      TH2   = TH2MN + (TH2MX-TH2MN)*Z(10)
      Q22   = AMGZ*TAN(TH2) + AMZ2
      Q2    = SQRT(Q22)
      WAT   = WAT*(TH2MX-TH2MN)*( (Q22-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Z1Z2.
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
      Q32MN = (AMZ+AMZ)**2
      Q32MX = (RS-AMZ)**2
      Q22   = AMZ2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(11)
      Q42   = AMZ2
      Q52   = AMZ2
C--
      WAT   = WATS*(Q32MX-Q32MN)
     .            *(xPI*AMZ*GMZ)**3
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
      CALL SGZZZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
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
      AMZZ2 = XPHASE(1,2)
      AMZZ  = SQRT(AMZZ2)
      EBM   = ROOTS/2
      EZ1   = QV(0,4)
      EZ2   = QV(0,5)
      EZ    = QV(0,2)
      RDDG  = 180/xPI
      COSZ  = XPHASE(2,1)
      PHIZ  = XPHASE(3,1)
      COSZ1 = XPHASE(2,2)
      PHIZ1 = XPHASE(3,2)
      CALL XHFILL( 'h00', DBLE(RS/ROOTS)    ,FUNC )
      CALL XHFILL( 'h01', DBLE(COSZ)        ,FUNC )
      CALL XHFILL( 'h02', DBLE(PHIZ*RDDG)   ,FUNC )
      CALL XHFILL( 'h03', DBLE(AMZZ/ROOTS)  ,FUNC )
      CALL XHFILL( 'h04', DBLE(COSZ1)       ,FUNC )
      CALL XHFILL( 'h05', DBLE(PHIZ1*RDDG)  ,FUNC )
      CALL XHFILL( 'h06', DBLE(Q2)          ,FUNC )
      CALL XHFILL( 'h07', DBLE(Q4)          ,FUNC )
      CALL XHFILL( 'h08', DBLE(Q5)          ,FUNC )
      CALL XHFILL( 'h09', DBLE(ICMB)        ,FUNC )
      CALL DHFILL( 'hd20', DBLE(EZ/EBM),  DBLE(COSZ)   , FUNC )
      CALL DHFILL( 'hd21', DBLE(EZ1/EBM), DBLE(EZ2/EBM), FUNC )
C--
C  That's it.
C--
      RETURN
      END
