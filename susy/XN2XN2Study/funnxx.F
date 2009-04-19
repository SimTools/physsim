C* (Update Record)
C*   92/06/09  K.Fujii	FUNC for e+ e- --> X0_2 X0_2 full amplitudes.
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
      PARAMETER  ( NP = 8, NHEL = 16 )
      INTEGER*4   IPV(3,2)
      REAL   *8   SG, PV(0:3,NP), AM(NP), AMX(2), GMX(2), AMR(3,3,2)
C--
      REAL   *8   AMSF(3,2,2), GTSF(3,2,2)
C--
      REAL   *8   QV(0:3,30)
      REAL   *8   WAT, WT
C--
      INTEGER*4   IHEL(NP,NHEL,2)
      DATA        (( IHEL(I,J,1), I=1,NP ), J=1,NHEL )  /
     .                   -1,+1,   -1,-1, 0,  -1,-1, 0,
     .                   -1,+1,   -1,-1, 0,  -1,+1, 0,
     .                   -1,+1,   -1,-1, 0,  +1,-1, 0,
     .                   -1,+1,   -1,-1, 0,  +1,+1, 0,
     .                   -1,+1,   -1,+1, 0,  -1,-1, 0,
     .                   -1,+1,   -1,+1, 0,  -1,+1, 0,
     .                   -1,+1,   -1,+1, 0,  +1,-1, 0,
     .                   -1,+1,   -1,+1, 0,  +1,+1, 0,
     .                   -1,+1,   +1,-1, 0,  -1,-1, 0,
     .                   -1,+1,   +1,-1, 0,  -1,+1, 0,
     .                   -1,+1,   +1,-1, 0,  +1,-1, 0,
     .                   -1,+1,   +1,-1, 0,  +1,+1, 0,
     .                   -1,+1,   +1,+1, 0,  -1,-1, 0,
     .                   -1,+1,   +1,+1, 0,  -1,+1, 0,
     .                   -1,+1,   +1,+1, 0,  +1,-1, 0,
     .                   -1,+1,   +1,+1, 0,  +1,+1, 0/
      DATA        (( IHEL(I,J,2), I=1,NP ), J=1,NHEL )  /
     .                   +1,-1,   -1,-1, 0,  -1,-1, 0,
     .                   +1,-1,   -1,-1, 0,  -1,+1, 0,
     .                   +1,-1,   -1,-1, 0,  +1,-1, 0,
     .                   +1,-1,   -1,-1, 0,  +1,+1, 0,
     .                   +1,-1,   -1,+1, 0,  -1,-1, 0,
     .                   +1,-1,   -1,+1, 0,  -1,+1, 0,
     .                   +1,-1,   -1,+1, 0,  +1,-1, 0,
     .                   +1,-1,   -1,+1, 0,  +1,+1, 0,
     .                   +1,-1,   +1,-1, 0,  -1,-1, 0,
     .                   +1,-1,   +1,-1, 0,  -1,+1, 0,
     .                   +1,-1,   +1,-1, 0,  +1,-1, 0,
     .                   +1,-1,   +1,-1, 0,  +1,+1, 0,
     .                   +1,-1,   +1,+1, 0,  -1,-1, 0,
     .                   +1,-1,   +1,+1, 0,  -1,+1, 0,
     .                   +1,-1,   +1,+1, 0,  +1,-1, 0,
     .                   +1,-1,   +1,+1, 0,  +1,+1, 0/
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
         IX     = 2
         JX     = 2
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
C      ( 4) : Z2 decay mode
C      ( 5) : e- helicity
C      ( 6) : m(ab)**2 for X0i
C      ( 7) : m(ac)**2
C      ( 8) : m(ab)**2 for X0j
C      ( 9) : m(ac)**2
C      (10) : m(X0i)**2
C      (11) : m(X0j)**2
C      (12) : cos(theta_X0i)
C      (13) : phi_X0i
C      (14) : cos(theta_a)     in X0i rest frame
C      (15) : phi_a            in X0i rest frame
C      (16) : phi_b            in X0i rest frame
C      (17) : cos(theta_a)     in X0j rest frame
C      (18) : phi_a            in X0j rest frame
C      (19) : phi_b            in X0j rest frame
C      (20) : E_beam spread of E-
C      (21) : E_beam spread of E+
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
c         CALL PROBEB(Z(20),Z(1),EBM,EMBM,EMBM0,PEM)
c         CALL PROBEB(Z(21),Z(2),EBM,EPBM,EPBM0,PEP)

         call jsfbeamgend(z(20), z(21), z(1), z(2), ebm, 
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
      ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(12)*1.D8
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
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
C-- 1st Z from 1st X0.
C
      XZ1 = BRVDK(IMD1MN-1,2) + (BRVDK(IMD1MX,2)-BRVDK(IMD1MN-1,2))*Z(4)
      DO 200 IMD = IMD1MN, IMD1MX
         IF ( XZ1.LT.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
      IMD = IMD1MX
C--
210   CONTINUE
      IDPTCL(1, 3) = MDVDK(1,IMD,2)
      IDPTCL(2, 3) = MDVDK(2,IMD,2)
      IDPTCL(3, 3) = MDVDK(3,IMD,2)
      IDPTCL(1, 4) = MDVDK(1,IMD,2)
      IDPTCL(2, 4) = MDVDK(2,IMD,2)
      IDPTCL(3, 4) = MDVDK(3,IMD,2)
      IDPTCL(1, 5) = 0
      IDPTCL(2, 5) = 0
      IDPTCL(3, 5) = 0
C--
      BRZ1 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD1 = IMD
      WAT  = WAT*(BRVDK(IMD1MX,2)-BRVDK(IMD1MN-1,2))/BRZ1
C-- 2nd Z from 2nd X0.
      XZ2  = (XZ1-BRVDK(IMD-1,2))/BRZ1
      XZ2 = BRVDK(IMD2MN-1,2) + (BRVDK(IMD2MX,2)-BRVDK(IMD2MN-1,2))*XZ2
      DO 300 IMD = IMD2MN, IMD2MX
         IF ( XZ2.LT.BRVDK(IMD,2) )              GO TO 310
300   CONTINUE
      IMD = 12
C--
310   CONTINUE
      IDPTCL(1, 6) = MDVDK(1,IMD,2)
      IDPTCL(2, 6) = MDVDK(2,IMD,2)
      IDPTCL(3, 6) = MDVDK(3,IMD,2)
      IDPTCL(1, 7) = MDVDK(1,IMD,2)
      IDPTCL(2, 7) = MDVDK(2,IMD,2)
      IDPTCL(3, 7) = MDVDK(3,IMD,2)
      IDPTCL(1, 8) = 0
      IDPTCL(2, 8) = 0
      IDPTCL(3, 8) = 0
C--
      BRZ2 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      IMD2 = IMD
      WAT  = WAT*(BRVDK(IMD2MX,2)-BRVDK(IMD2MN-1,2))/BRZ2
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
C     (1,2,3,4,5,6,7,8) = (e-,e+, f,fb,X0_1, f,fb,X0_1)
C--
      AM(1)   = AMF(IDPTCL(1,1),IDPTCL(2,1),IDPTCL(3,1))
      AM(2)   = AMF(IDPTCL(1,2),IDPTCL(2,2),IDPTCL(3,2))
      AM(3)   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM(4)   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM(5)   = SZM(1)
      AM(6)   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AM(7)   = AMF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AM(8)   = SZM(1)
C--
      AMR(1,1,1) = AMSF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AMR(2,1,1) = GTSF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AMR(3,1,1) = AMR(1,1,1) + AM(4)
      AMR(1,2,1) = AMSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AMR(2,2,1) = GTSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AMR(3,2,1) = AMR(1,2,1) + AM(3)
      AMR(1,3,1) = AMZ
      AMR(2,3,1) = GMZTOT
      AMR(3,3,1) = AMR(1,3,1) + AM(5)
C--
      AMR(1,1,2) = AMSF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AMR(2,1,2) = GTSF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AMR(3,1,2) = AMR(1,1,2) + AM(7)
      AMR(1,2,2) = AMSF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AMR(2,2,2) = GTSF(IDPTCL(1,7),IDPTCL(2,7),IDPTCL(3,7))
      AMR(3,2,2) = AMR(1,2,2) + AM(6)
      AMR(1,3,2) = AMZ
      AMR(2,3,2) = GMZTOT
      AMR(3,3,2) = AMR(1,3,2) + AM(8)
C--
C  Set phase space variables.
C--
      CALL GETXPH(RS,AM,AMX,GMX,AMR,Z(6),IPV,XPHASE,WT)
      IF ( WT.LE.0.D0 ) THEN
         FUNC = 0
         RETURN
      ENDIF
      WAT = WAT*WT
C--
C  Calculate differential cross sections.
C--
      SG = 0
      DO 400 IH1 = 1, 2
         IHLCMB(5) = (-1)**IH1 
         DO 40 IH2 = 1, 2
            IHLCMB(8) = (-1)**IH2
            CALL SGNXXF(IDPTCL,IHLCMB,AM,IPV,XPHASE,SGH,PV)
            SG = SG + SGH
40       CONTINUE
400   CONTINUE
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
      CALL PSUMxx(PV(0,5),PV(0,8),QV(0,1))
      PTOP  = SQRT( QV(1,1)**2 + QV(2,1)**2 )
C--      
      COSL1  =  PV(3,3)/UABV3D(PV(1,3))
      COSLB1 = -PV(3,4)/UABV3D(PV(1,4))
      COSL2  =  PV(3,6)/UABV3D(PV(1,6))
      COSLB2 = -PV(3,7)/UABV3D(PV(1,7))
      COST   = XPHASE(2,1)
      PHIT   = XPHASE(3,1)
      COSBB1 = XPHASE(2,2)
      PHIBB1 = XPHASE(3,2)
      COSBB2 = XPHASE(2,4)
      PHIBB2 = XPHASE(3,4)
      COMB   = JCMB + NHEL*(ICMB-1) + .1D0
      QXI2   = XPHASE(1,2)
      Q32    = XPHASE(1,3)
      Q42    = XPHASE(2,3)
      QXJ2   = XPHASE(1,4)
      Q62    = XPHASE(1,5)
      Q72    = XPHASE(2,5)
C--      
      CALL XHFILL( 'h01', COST               ,FUNC )
      CALL XHFILL( 'h02', PHIT*xRD2DG        ,FUNC )
      CALL XHFILL( 'h03', COSBB1             ,FUNC )
      CALL XHFILL( 'h04', PHIBB1*xRD2DG      ,FUNC )
      CALL XHFILL( 'h05', SQRT(Q32)          ,FUNC )
      CALL XHFILL( 'h06', SQRT(Q42)          ,FUNC )
      CALL XHFILL( 'h07', COSBB2             ,FUNC )
      CALL XHFILL( 'h08', PHIBB2*xRD2DG      ,FUNC )
      CALL XHFILL( 'h09', SQRT(Q62)          ,FUNC )
      CALL XHFILL( 'h10', SQRT(Q72)          ,FUNC )
      CALL XHFILL( 'h11', PTOP               ,FUNC )
      CALL XHFILL( 'h12', COMB               ,FUNC )
      CALL XHFILL( 'h13', PV(0,3)            ,FUNC )
      CALL XHFILL( 'h14', PV(0,4)            ,FUNC )
      CALL XHFILL( 'h15', PV(0,6)            ,FUNC )
      CALL XHFILL( 'h16', PV(0,7)            ,FUNC )
      CALL XHFILL( 'h17', COSL1              ,FUNC )
      CALL XHFILL( 'h18', COSLB1             ,FUNC )
      CALL XHFILL( 'h19', COSL2              ,FUNC )
      CALL XHFILL( 'h20', COSLB2             ,FUNC )
      CALL XHFILL( 'h21', RS/ROOTS           ,FUNC )
      CALL XHFILL( 'h22', SQRT(QXI2)         ,FUNC )
      CALL XHFILL( 'h23', SQRT(QXJ2)         ,FUNC )
      CALL XHFILL( 'h24', IMD1+.1D0          ,FUNC )
      CALL XHFILL( 'h24', IMD2+12.1D0        ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
