C* (Update Record)
C*   2010/06/20  J.Tian & K.Fujii	FUNC for e+ e- --> nu nubar Z H full amplitudes.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "NNZHCNF.inc"
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
#include "USRPRM.inc"
      PARAMETER  ( NP = 7, NHEL = 2 )
      REAL   *8   SG, PV(0:3,NP)
C--
      REAL   *8   QV(0:3,5)
      REAL   *8   Q3, Q32, Q3MN, Q3MX
      REAL   *8   Q4, Q42, Q4MN, Q4MX, TH4MN, TH4MX, TH4
      REAL   *8   WAT
C--
      INTEGER*4  IHEL(NP,NHEL,2)
      DATA ( ( IHEL(I,J,1),I=1,NP ), J= 1,NHEL )  /
     .             -1,+1,   -1, +1,   +1,-1,0,
     .             -1,+1,   -1, +1,   -1,+1,0/
      DATA ( ( IHEL(I,J,2),I=1,NP ), J= 1,NHEL )  /
     .             +1,-1,   -1, +1,   +1,-1,0,
     .             +1,-1,   -1, +1,   -1,+1,0/
C--
#ifdef WITH_DBD_STANDARD
      REAL*8 EPSISR,FISR3,FISR4,Z3INI,Z4INI 
      REAL*8 PISRDBL(0:3,2)
#endif
#ifdef DEBUG      
      INTEGER*4  JDP(2,NP)
#endif
#if 1
      common /event_number/ nseqevent
#endif
      DATA NCALL /  0 /
      EXTERNAL   XHFILL, DHFILL
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
C--
         AMW2  = AMW*AMW
         GMZ   = GMZTOT
         AMGZ  = AMZ*GMZ
         AMZ2  = AMZ*AMZ
         AMGZ2 = AMGZ*AMGZ
C--
         AME   = AMF(1,2,1)
C--
         xALF0 = 1/137.0359895E0
C--
         XZ1MN = BRVDK(IMDMN-1,2)
         XZ1MX = BRVDK(IMDMX  ,2)
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
C      ( 3) : ISR for e-
C      ( 4) : ISR for e+ ( or e-/e+ selection when ISR_LLR_Order=0)
C      ( 5) : helicity
C      ( 6) : xi   (E_nu)
C      ( 7) : m(Z)**2
C      ( 8) : eta  (cos_nu)
C      ( 9) : zeta (cos_nb) in (nb Z H) frame
C      (10) : m(ZH)**2
C      (11) : cos_Z        in (Z H)  frame
C      (12) : phi_Z        in (Z H)  frame
C      (13) : phi_nu
C      (14) : phi_nb  	   in (nb Z H) frame
C      (15) : cos_fb       in Z      frame
C      (16) : phi_f        in Z      frame
C      (17) : final state combination.
C      (18) : E_beam spread of E-
C      (19) : E_beam spread of E+
C      (20) : random variable for Pt of ISR photon from e- beam
C      (21) : random variable for Pt of ISR photon from e+ beam
C      (22) : random variable for azimuthal angle of ISR photon from e- beam
C      (23) : random variable for azimuthal angle of ISR photon from e+ beam
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Beam energy.
C--
      EBM       = ROOTS/2
      RS        = ROOTS
      EMBM      = RS/2
      EPBM      = RS/2
      QED       = 1
C--
C  First decide e- and e+ energy.
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
         EBM = ROOTS/2
         call jsfbeamgenr(Z(1), Z(2), Z(18), Z(19), ebm,
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
#else 
C -- when WITH_DBD_STANDARD is set
      IF ( ISRBM.GE.3 ) THEN
        IF ( ISRBM.EQ.3 ) THEN
          call jsfbeamgenr(Z(1), Z(2), Z(18), Z(19), ebm,
     >                    embm, epbm, embm0, epbm0, dpdebm )
        ELSE
          call whizard_spectrum( isrbm, ROOTS, z(1), z(2),
     >       ebm, dpdebm, embm, epbm )
          embm0=embm
          epbm0=epbm
        ENDIF
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
#endif
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
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
#ifdef WITH_DBD_STANDARD
          CALL ISR_WHIZARD_SPECTRUM(Z(3),Z(4),Z(20),Z(21),Z(22),Z(23),
     >      xALF0, ISR_LLA_ORDER, ROOTS, EMBM, EPBM, QV(0,1),
     >     PISRDBL, WAT, RS)
           IF( RS < 0.0 ) THEN
              FUNC=0
              RETURN
           ENDIF
           DO 10500 KK=0,3
              PVOUT(KK,10)=PISRDBL(KK,1)
              PVOUT(KK,11)=PISRDBL(KK,2)
10500      CONTINUE
#endif
        ENDIF
      ENDIF
C      print *,"Debug_W: ",ROOTS,EMBM,EPBM,PVOUT(0,10)
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
      IDPTCL(1, 3) = 1
      IDPTCL(2, 3) = 1
      IDPTCL(3, 3) = 1
      IDPTCL(1, 4) = 1
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = 1
C-- Z.
      XZ1 = XZ1MN + Z(17)*(XZ1MX-XZ1MN)
C      print *,"Debug_IMD: ",XZ1MN,XZ1MX,XZ1,Z(17),IMDMN,IMDMX
      DO 200 IMD = IMDMN, IMDMX
         IF ( XZ1.LT.BRVDK(IMD,2) )              GO TO 210
200   CONTINUE
C--
210   IDPTCL(1,5) = MDVDK(1,IMD,2)
      IDPTCL(2,5) = MDVDK(2,IMD,2)
      IDPTCL(3,5) = MDVDK(3,IMD,2)
      IDPTCL(1,6) = MDVDK(1,IMD,2)
      IDPTCL(2,6) = MDVDK(2,IMD,2)
      IDPTCL(3,6) = MDVDK(3,IMD,2)
C--
      BRZ1 = BRVDK(IMD,2) - BRVDK(IMD-1,2)
      WAT  = WAT*(XZ1MX-XZ1MN)/BRZ1
C--
C  Select helicity combination.
C--
      XHL  = Z(5)
      HLM  = (1-POLE)/2
      IF ( XHL.LT.HLM ) THEN
         ICMB = 1
         JCMB = NHEL*XHL/HLM + 1
      ELSE
         ICMB = 2
         JCMB = NHEL*(XHL-HLM)/(1-HLM) + 1
      ENDIF
      WAT  = WAT*NHEL
      JCMB = MIN(JCMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,JCMB,ICMB),IHLCMB(1))
      XHCMB = NHEL*(ICMB-1) + JCMB
C--
C  Set kinematics.
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMF(IDPTCL(1,5),IDPTCL(2,5),IDPTCL(3,5))
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
      AM7   = AMH
C      print *,"Debug_M: ",AM3,AM4,AM5,AM6,AM7
C-- q^2_ZH.
      Q32MN = (AM5 + AM6 + AM7)**2
      Q32MX = (RS - ( AM3 + AM4 ))**2
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(10)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
C-- q^2_Z.
      Q4MN  = AM5 + AM6
      Q4MX  = Q3 - AM7
      TH4MN = ATAN((Q4MN-AMZ)*(Q4MN+AMZ)/AMGZ)
      TH4MX = ATAN((Q4MX-AMZ)*(Q4MX+AMZ)/AMGZ)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(7)
      Q42   = AMGZ*TAN(TH4) + AMZ2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMZ2)**2 + AMGZ2 )/AMGZ
C-- Xi_lo.
      XILO  = LOG( (Q3-AM3)*(Q3+AM3)/S )
      XIHI  = LOG( 1 - 2*AM3/RS )
      XI    = XILO + (XIHI-XILO)*Z(6)
      WAT   = WAT*(XIHI-XILO)
      EXPXI = EXP(XI)
C-- Eta.
      DM1   = AMW2/S
      DP1   = 1
      ETLO  = -LOG( (1+DM1)/DP1 )/2
      ETHI  =  LOG( (1+DP1)/DM1 )/2
      ETA   = ETLO + (ETHI-ETLO)*Z(8)
      WAT   = WAT*(ETHI-ETLO)
C      print *,"Debug_P: ",Q3,Q4,XI,ETA
C--
C  Set XPHASE.
C--
      XPHASE(1,1) = S
      XPHASE(2,1) = ETA
      XPHASE(3,1) = Z(13)*x2PI
      XPHASE(1,2) = XI
      XPHASE(2,2) = Z(9)
      XPHASE(3,2) = Z(14)*x2PI
      XPHASE(1,3) = Q32
      XPHASE(2,3) = -1 + 2*Z(11)
      XPHASE(3,3) = Z(12)*x2PI
      XPHASE(1,4) = Q42
      XPHASE(2,4) = -1 + 2*Z(15)
      XPHASE(3,4) = Z(16)*x2PI
C--
      WAT   = WAT*4*(x2PI)**4
C--
C  Calculate differential cross sections.
C--
      CALL SGNNZHF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG     = SG*QED
      FUNC   = SG*WAT
C      print *,"Debug_SG: ",SG,FUNC
      IF ( FUNC.LE.0.D0 ) RETURN
C--
C  Boost everything to laboratroy frame.
C--
#if 0
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
#endif
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
      RDDG   = 180/xPI
      QWN2   = S*EXPXI + AM3*AM3
      QWN    = SQRT(QWN2)
      EXP2ET = EXP(2*XPHASE(2,1))
      SH1    = MIN( SQRT( (1+DP1+DM1)/(1+EXP2ET) - DM1 ), 1.D0 )
      DLT    = 2*AMZ*AMZ/((QWN-AMZ)*(QWN+AMZ))
      OPDZT  = ( 1 + 1/DLT )**XPHASE(2,2)
      SH2    = MIN( SQRT( DLT*( OPDZT - 1 ) ), 1.D0 )
C--
      XE    = 1 + AM3*AM3/S - QWN2/S
      COSE  = ( 1 - xSQRT2*SH1 )*( 1 + xSQRT2*SH1 )
      PHIE  = XPHASE(3,1)
      COSN  = ( 1 - xSQRT2*SH2 )*( 1 + xSQRT2*SH2 )
      PHIN  = XPHASE(3,2)
      QW    = Q3
      CSDA  = XPHASE(2,3)
      PHID  = XPHASE(3,3)
C--
      QV(0,2) = PV(0,3) - PV(0,1)
      QV(1,2) = PV(1,3) - PV(1,1)
      QV(2,2) = PV(2,3) - PV(2,1)
      QV(3,2) = PV(3,3) - PV(3,1)
      AQV     = SQRT(QV(1,2)**2+QV(2,2)**2+QV(3,2)**2)
      TEE     = (QV(0,2)-AQV)*(QV(0,2)+AQV)
C--
      QV(0,5) = PV(0,5) + PV(0,6)
      QV(1,5) = PV(1,5) + PV(1,6)
      QV(2,5) = PV(2,5) + PV(2,6)
      QV(3,5) = PV(3,5) + PV(3,6)
      XTW     = 2*SQRT(QV(1,5)**2+QV(2,5)**2+QV(3,5)**2)/ROOTS
C--
      XEL     = 2*PV(0,3)/ROOTS
      XNU     = 2*PV(0,4)/ROOTS
      XWP     = 2*QV(0,5)/ROOTS
C--
      CALL XHFILL( 'h01', COSE              ,FUNC )
      CALL XHFILL( 'h02', PHIE*RDDG         ,FUNC )
      CALL XHFILL( 'h03', XE                ,FUNC )
      CALL XHFILL( 'h04', Q4                ,FUNC )
      CALL XHFILL( 'h05', COSN              ,FUNC )
      CALL XHFILL( 'h06', PHIN*RDDG         ,FUNC )
      CALL XHFILL( 'h07', CSDA              ,FUNC )
      CALL XHFILL( 'h08', PHID*RDDG         ,FUNC )
      CALL XHFILL( 'h09', RS/ROOTS          ,FUNC )
      CALL XHFILL( 'h10', IMD+.1D0          ,FUNC )
      CALL XHFILL( 'h11', XHCMB+.1D0        ,FUNC )
      CALL XHFILL( 'h12', ETA               ,FUNC )
      CALL XHFILL( 'h13', XTW               ,FUNC )
      CALL XHFILL( 'h14', XEL               ,FUNC )
      CALL XHFILL( 'h15', XNU               ,FUNC )
      CALL XHFILL( 'h16', XWP               ,FUNC )
      CALL XHFILL( 'h17', QW                ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
