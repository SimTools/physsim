C* (Update Record)
C*  95/03/12  K.Fujii		Converted to unix. 
C*				To avoid problems related to kinematic
C*				boundaries, we introduced a tiny offset
C*				EPSL. Now use INSMCN instead of CONINT.

      REAL FUNCTION FUNC*8(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "ZZCNSF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6, NHEL = 4 )
      REAL   *4   SG, PV(0:3,NP), QV(0:3,0:5)
      REAL   *8   WAT
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
C
C                              Z        Z
C                            _        _
C                  e- e+     f   f    f   f
C               ----------------------------------
     .             -1,+1,   -1, +1,  -1, +1,
     .             -1,+1,   -1, +1,  +1, -1,
     .             -1,+1,   +1, -1,  -1, +1,
     .             -1,+1,   +1, -1,  +1, -1,
     .             +1,-1,   -1, +1,  -1, +1,
     .             +1,-1,   -1, +1,  +1, -1,
     .             +1,-1,   +1, -1,  -1, +1,
     .             +1,-1,   +1, -1,  +1, -1 /
      DATA NCALL / 0 /
C--
      DATA EPSL / 1.0 /
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
C--
         GMZ   = GMZTOT
         AMGZ  = AMZ*GMZ
         AMZ2  = AMZ*AMZ
         AMGZ2 = AMGZ*AMGZ
C--
         AME   = AMF(1,2,1)
         xALF0 = 1/137.0359895E0
C--
         XZ1MN = BRVDK(IMD1MN-1,2)
         XZ1MX = BRVDK(IMD1MX  ,2)
         XZ2MN = BRVDK(IMD2MN-1,2)
         XZ2MX = BRVDK(IMD2MX  ,2)
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
C      ( 2) : e- beam gaussian spread
C      ( 3) : e+ beam
C      ( 4) : e+ beam gaussian spread
C      ( 5) : bremsstrahlng
C      ( 6) : e- helicity
C      ( 7) : helicity combination for final states.
C      ( 8) : m(Z)**2
C      ( 9) : m(Z)**2
C      (10) : cos(theta_Z)
C      (11) : phi_Z
C      (12) : cos(theta_fb)     in Z rest frame
C      (13) : phi_fb            in Z rest frame
C      (14) : cos(theta_fb)     in Z rest frame
C      (15) : phi_fb            in Z rest frame
C      (16) : final state combination
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
         ZG  = Z(5)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(2)*1.E4
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
C--
      XZ1 = XZ1MN + Z(16)*(XZ1MX-XZ1MN)
      DO 100 IMD = IMD1MN, IMD1MX
         IF ( XZ1.LT.BRVDK(IMD,2) )              GO TO 110
100   CONTINUE
C--
110   IDPTCL(1, 3) = MDVDK(1,IMD,2)
      IDPTCL(2, 3) = MDVDK(2,IMD,2)
      IDPTCL(3, 3) = MDVDK(3,IMD,2)
      IDPTCL(1, 4) = MDVDK(1,IMD,2)
      IDPTCL(2, 4) = MDVDK(2,IMD,2)
      IDPTCL(3, 4) = MDVDK(3,IMD,2)
C--
      BRZ1 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD1 = IMD
      WAT  = WAT*(XZ1MX-XZ1MN)/BRZ1
C--
      XZ2  = XZ2MN + ((XZ1-BRVDK(IMD-1,2))/BRZ1)*(XZ2MX-XZ2MN)
      DO 200 IMD = IMD2MN, IMD2MX
         IF ( XZ2.LT.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1, 5) = MDVDK(1,IMD,2)
      IDPTCL(2, 5) = MDVDK(2,IMD,2)
      IDPTCL(3, 5) = MDVDK(3,IMD,2)
      IDPTCL(1, 6) = MDVDK(1,IMD,2)
      IDPTCL(2, 6) = MDVDK(2,IMD,2)
      IDPTCL(3, 6) = MDVDK(3,IMD,2)
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
      ICMB = JHL + IOFF
      IF ( IDPTCL(2,3).NE.1 .OR. IDPTCL(3,3).NE.1 ) THEN
         IF ( IDPTCL(2,5).NE.1 .OR. IDPTCL(3,5).NE.1 ) THEN
            WAT  = WAT*NHEL
         ELSE
            ICMB = 2*((ICMB+1)/2)
            WAT  = WAT*NHEL/2
         ENDIF
      ELSE
         IF ( IDPTCL(2,5).NE.1 .OR. IDPTCL(3,5).NE.1 ) THEN
            ICMB = 4*((ICMB-1)/4) + MOD(ICMB,2) + 3
            WAT  = WAT*NHEL/2
         ELSE
            ICMB = 4*((ICMB-1)/4) + 4
            WAT  = WAT*NHEL/4
         ENDIF
      ENDIF
      CALL UVCOPY(NP,IHEL(1,ICMB),IHLCMB)
C--
C  Set invariant masses squared.
C--
      AM3   = AMF(IDPTCL(1, 3),IDPTCL(2, 3),IDPTCL(3, 3))
      AM4   = AMF(IDPTCL(1, 4),IDPTCL(2, 4),IDPTCL(3, 4))
      AM5   = AMF(IDPTCL(1, 5),IDPTCL(2, 5),IDPTCL(3, 5))
      AM6   = AMF(IDPTCL(1, 6),IDPTCL(2, 6),IDPTCL(3, 6))
C-- Z.
      Q2MN  = AM3 + AM4 + EPSL
      Q2MX  = MAX( RS - (AM5+AM6) - EPSL, Q2MN )
      TH2MN = ATAN((Q2MN-AMZ)*(Q2MN+AMZ)/AMGZ)
      TH2MX = ATAN((Q2MX-AMZ)*(Q2MX+AMZ)/AMGZ)
      TH2   = TH2MN + (TH2MX-TH2MN)*Z(8)
      Q22   = AMGZ*TAN(TH2) + AMZ2
      Q2    = SQRT(Q22)
      WAT   = WAT*(TH2MX-TH2MN)*( (Q22-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Z.
      Q3MN  = AM5 + AM6 + EPSL
      Q3MX  = MAX( RS - Q2 - EPSL, Q3MN )
      TH3MN = ATAN((Q3MN-AMZ)*(Q3MN+AMZ)/AMGZ)
      TH3MX = ATAN((Q3MX-AMZ)*(Q3MX+AMZ)/AMGZ)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(9)
      Q32   = AMGZ*TAN(TH3) + AMZ2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMZ2)**2 + AMGZ2 )/AMGZ
C>>>
CDBG        PRINT *, ' AMZ, GMZ, Q3 = ', AMZ, GMZ, Q3
C>>>
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q22
      XPHASE(1,3) = Q32
C--
C  Angular variables.
C--
      DO 400 IB = 1, 3
         XPHASE(2,IB) = -1 + 2*Z(2*IB+8)
         XPHASE(3,IB) = x2PI  *Z(2*IB+9)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Identical particle factor.
C--
      IDNT = 2
      WAT  = WAT/IDNT
C--
C  Calculate differential cross sections.
C--
      CALL SGZZF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG   = SG*QED
      FUNC = SG*WAT
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
CDBG       PRINT *, ' QV = ', (QV(K,1),K=0,3)
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
      CALL PSUMxx(PV(0,5),PV(0,6),QV(0,3))
      CALL UVZERO(4,QV(0,0))
      QV(0,0) = ROOTS
      CALL PSUBxx(QV(0,0),QV(0,3),QV(0,4))
      CALL RDOTxx(QV(0,4),QV(0,4),QMIS2)
      IF ( QMIS2.GE.0. ) THEN
         QMIS  = SQRT(QMIS2)
      ELSE
         QMIS  = -SQRT(-QMIS2)
      ENDIF
C--
      RDDG  = 180/xPI
      COSZ  = XPHASE(2,1)
      PHIZ  = XPHASE(3,1)
      COSM  = XPHASE(2,2)
      PHIM  = XPHASE(3,2)
      COSF  = XPHASE(2,3)
      PHIF  = XPHASE(3,3)
      CALL XHFILL(  1, DBLE(RS  /ROOTS)  ,FUNC )
      CALL XHFILL(  2, DBLE(Q2)          ,FUNC )
      CALL XHFILL(  3, DBLE(Q3)          ,FUNC )
      CALL XHFILL(  4, DBLE(QMIS/ROOTS)  ,FUNC )
      CALL XHFILL(  5, DBLE(COSZ)        ,FUNC )
      CALL XHFILL(  6, DBLE(PHIZ)*RDDG   ,FUNC )
      CALL XHFILL(  7, DBLE(COSF)        ,FUNC )
      CALL XHFILL(  8, DBLE(PHIF)*RDDG   ,FUNC )
      CALL XHFILL(  9, DBLE(COSM)        ,FUNC )
      CALL XHFILL( 10, DBLE(PHIM)*RDDG   ,FUNC )
      CALL XHFILL( 11, DBLE(ICMB)+.1     ,FUNC )
      CALL XHFILL( 12, DBLE(IMD1)+.1     ,FUNC )
      CALL XHFILL( 12, DBLE(IMD2)+12.1   ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
