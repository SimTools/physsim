C* (Update Record)
C*    92/06/09  K.Fujii	  FUNC for e+ e- --> X0_1 X0_2 full amplitudes.
C*    03/09/10  K.Fujii	  Modified to allow
C*                        Z2->snL+nu, seR+e, suR+u, sdR+d (1st and 2nd)
C*                          ->snL+nu, se1+e, su1+u, sd1+d (3rd gen)
C*                        !!! CAUTION !!!
C*                          seL, se2, suL, sdL decays are ignored.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "XNXNCN.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6, NHEL = 16 )
      INTEGER*4   IPV(3,2)
      REAL   *8   SG, PV(0:3,NP), AM(NP), AMX(2), GMX(2), AMR(3,3)
C--
      REAL   *8   AMSF(3,2,2), GTSF(3,2,2)
C--
      REAL   *8   QV(0:3,30)
      REAL   *8   WAT, WT
C--
      INTEGER*4   IHEL(NP,NHEL,2)
      DATA        (( IHEL(I,J,1), I=1,NP ), J=1,NHEL )  /
     .                   -1,+1,   -1,  -1, -1,-1,
     .                   -1,+1,   -1,  -1, -1,+1,
     .                   -1,+1,   -1,  -1, +1,-1,
     .                   -1,+1,   -1,  -1, +1,+1,
     .                   -1,+1,   -1,  +1, -1,-1,
     .                   -1,+1,   -1,  +1, -1,+1,
     .                   -1,+1,   -1,  +1, +1,-1,
     .                   -1,+1,   -1,  +1, +1,+1,
     .                   -1,+1,   +1,  -1, -1,-1,
     .                   -1,+1,   +1,  -1, -1,+1,
     .                   -1,+1,   +1,  -1, +1,-1,
     .                   -1,+1,   +1,  -1, +1,+1,
     .                   -1,+1,   +1,  +1, -1,-1,
     .                   -1,+1,   +1,  +1, -1,+1,
     .                   -1,+1,   +1,  +1, +1,-1,
     .                   -1,+1,   +1,  +1, +1,+1/
      DATA        (( IHEL(I,J,2), I=1,NP ), J=1,NHEL )  /
     .                   +1,-1,   -1,  -1, -1,-1,
     .                   +1,-1,   -1,  -1, -1,+1,
     .                   +1,-1,   -1,  -1, +1,-1,
     .                   +1,-1,   -1,  -1, +1,+1,
     .                   +1,-1,   -1,  +1, -1,-1,
     .                   +1,-1,   -1,  +1, -1,+1,
     .                   +1,-1,   -1,  +1, +1,-1,
     .                   +1,-1,   -1,  +1, +1,+1,
     .                   +1,-1,   +1,  -1, -1,-1,
     .                   +1,-1,   +1,  -1, -1,+1,
     .                   +1,-1,   +1,  -1, +1,-1,
     .                   +1,-1,   +1,  -1, +1,+1,
     .                   +1,-1,   +1,  +1, -1,-1,
     .                   +1,-1,   +1,  +1, -1,+1,
     .                   +1,-1,   +1,  +1, +1,-1,
     .                   +1,-1,   +1,  +1, +1,+1/
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
         IX     = 1
         JX     = 2
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
C      ( 4) : e- helicity
C      ( 5) : Z2 decay mode
C      ( 6) : m(ab)**2
C      ( 7) : m(ac)**2
C      ( 8) : m(X0j)**2
C      ( 9) : cos(theta_X0j)
C      (10) : phi_X0j
C      (11) : cos(theta_a)     in X0j rest frame
C      (12) : phi_a            in X0j rest frame
C      (13) : phi_b            in X0j rest frame
C      (14) : E_beam spread of E-
C      (15) : E_beam spread of E+
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
          CALL PROBEB(Z(14),Z(1),EBM,EMBM,EMBM0,PEM)
         CALL PROBEB(Z(15),Z(2),EBM,EPBM,EPBM0,PEP)
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
         DPDEBM    = PEM*PEP
      ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(9)*1.D4
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
      IDPTCL(1, 3) = 0
      IDPTCL(2, 3) = 0
      IDPTCL(3, 3) = 0
      IDPTCL(1, 6) = 0
      IDPTCL(2, 6) = 0
      IDPTCL(3, 6) = 0
C-- Z from X0.
      XZ = Z(5)
      DO 200 IMD = 1, 12
         IF ( XZ.LT.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1, 4) = MDVDK(1,IMD,2)
      IDPTCL(2, 4) = MDVDK(2,IMD,2)
      IDPTCL(3, 4) = MDVDK(3,IMD,2)
      IDPTCL(1, 5) = MDVDK(1,IMD,2)
      IDPTCL(2, 5) = MDVDK(2,IMD,2)
      IDPTCL(3, 5) = MDVDK(3,IMD,2)
C--
      BRZ = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      WAT  = WAT/BRZ
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( Z(4).LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*Z(4)/HLM + 1
      ELSE
         ICMB = 2
         JCMB = NHEL*(Z(4)-HLM)/(1-HLM) + 1
      ENDIF
      WAT = WAT*NHEL
      JCMB = MIN(JCMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
C--
C  Select kinematic variables.
C     (1,2,3,4,5,6) = (e-,e+, X0_1, f,fb,X0_1)
C--
      AM(1)   = AMF(IDPTCL(1,1),IDPTCL(2,1),IDPTCL(3,1))
      AM(2)   = AMF(IDPTCL(1,2),IDPTCL(2,2),IDPTCL(3,2))
      AM(3)   = SZM(1)
      AM(4)   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM(5)   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM(6)   = SZM(1)
C--
      AMR(1,1) = AMSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4)) 
      AMR(2,1) = GTSF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AMR(3,1) = AMR(1,1) + AM(5)
      
      AMR(1,2) = AMSF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AMR(2,2) = GTSF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AMR(3,2) = AMR(1,2) + AM(4)
      AMR(1,3) = AMZ
      AMR(2,3) = GMZTOT
      AMR(3,3) = AMR(1,3) + AM(6)
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
      CALL SGNXXF(IDPTCL,IHLCMB,AM,IPV,XPHASE,SG,PV)
      SG     = SG*QED
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
      CALL PSUMxx(PV(0,3),PV(0,6),QV(0,1))
      PTOP  = SQRT( QV(1,1)**2 + QV(2,1)**2 )
C--
      CALL PSUMxx(PV(0,4),PV(0,5),QV(0,2))
      CALL PSUMxx(QV(0,2),PV(0,6),QV(0,12))
      CALL PMIRRx(QV(0,12),QV(0,12))
      CALL BOOSTx(PV(0,4),QV(0,12),QV(0,14))
      CALL BOOSTx(PV(0,5),QV(0,12),QV(0,15))
      CS14 = UDOT3D(QV(1,14),QV(1,12))/UABV3D(QV(1,14))/UABV3D(QV(1,12))
      CS15 = UDOT3D(QV(1,15),QV(1,12))/UABV3D(QV(1,15))/UABV3D(QV(1,12))
C--
      CALL PSUMxx(PV(0,6),PV(0,4),QV(0,21))
      CALL PMIRRx(QV(0,21),QV(0,21))
      CALL BOOSTx(PV(0,4),QV(0,21),QV(0,24))
      CS24 = UDOT3D(QV(1,24),QV(1,21))/UABV3D(QV(1,24))/UABV3D(QV(1,21))
      CALL PSUMxx(PV(0,6),PV(0,5),QV(0,22))
      CALL PMIRRx(QV(0,22),QV(0,22))
      CALL BOOSTx(PV(0,5),QV(0,22),QV(0,25))
      CS25 = UDOT3D(QV(1,25),QV(1,22))/UABV3D(QV(1,25))/UABV3D(QV(1,22))
C--      
      COSL  =  PV(3,4)/UABV3D(PV(1,4))
      COSLB = -PV(3,5)/UABV3D(PV(1,5))
      PT45  = MAX(SQRT( QV(1,2)**2 + QV(2,2)**2 ),.1D0)
      PT4   = MAX(SQRT( PV(1,4)**2 + PV(2,4)**2 ),.1D0)
      PT5   = MAX(SQRT( PV(1,5)**2 + PV(2,5)**2 ),.1D0)
      ACOP  = ( PV(1,4)*PV(1,5) + PV(2,4)*QV(2,5) )/PT4/PT5
      IF ( ABS(ACOP).GE.1.0D0 ) ACOP = SIGN(1.0,ACOP)
      ACOP  = xPI - ACOS(ACOP)
      COST  = XPHASE(2,1)
      PHIT  = XPHASE(3,1)
      COSBB = XPHASE(2,2)
      PHIBB = XPHASE(3,2)
      COMB  = JCMB + NHEL*(ICMB-1) + .1D0
      QXJ2  = XPHASE(1,2)
      Q42   = XPHASE(1,3)
      Q52   = XPHASE(2,3)
C--      
      CALL XHFILL( 'h01', COST              ,FUNC )
      CALL XHFILL( 'h02', PHIT*xRD2DG       ,FUNC )
      CALL XHFILL( 'h03', COSBB             ,FUNC )
      CALL XHFILL( 'h04', PHIBB*xRD2DG      ,FUNC )
      CALL XHFILL( 'h05', SQRT(Q42)         ,FUNC )
      CALL XHFILL( 'h06', SQRT(Q52)         ,FUNC )
      CALL XHFILL( 'h07', PTOP              ,FUNC )
      CALL XHFILL( 'h08', ACOP*xRD2DG       ,FUNC )
      CALL XHFILL( 'h09', COMB              ,FUNC )
      CALL XHFILL( 'h10', PV(0,4)           ,FUNC )
      CALL XHFILL( 'h11', PV(0,5)           ,FUNC )
      CALL XHFILL( 'h12', COSL              ,FUNC )
      CALL XHFILL( 'h13', COSLB             ,FUNC )
      CALL XHFILL( 'h14', RS/ROOTS          ,FUNC )
      CALL XHFILL( 'h15', CS14              ,FUNC )
      CALL XHFILL( 'h16', CS15              ,FUNC )
      CALL XHFILL( 'h17', SQRT(QXJ2)        ,FUNC )
      CALL XHFILL( 'h18', CS24              ,FUNC )
      CALL XHFILL( 'h19', CS25              ,FUNC )
      CALL XHFILL( 'h20', DBLE(IMD)+0.1D0   ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
