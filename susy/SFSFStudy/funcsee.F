C*  92/06/29  K.Fujii	FUNC for e+ e- --> SE+ SE- full amplitudes.
C*  95/03/21  K.Fujii	QED vtx correction ( 1+3*BTE/4 ) and
C*                      ( 1 + 2*ALF0/PI * ( PI**2/6 - 1/4 ) )
C*                      included.
C*  95/05/05  K.Fujii	Random number shuffling and ISR and BMEFF
C*			treatment on ISRBM.
C*  00/06/24  K.Fujii	Modified for Bases++.
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SFSFCN.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
      PARAMETER  ( NP = 6 )
      REAL   *4   SG, PV(0:3,NP)
C--
      REAL   *4   QV(0:3,10)
C--
      PARAMETER   ( NHEL = 4 )
      INTEGER*4   IHEL(NP,NHEL,3,2)
      DATA        ((( IHEL(I,J,K,1),I=1,NP),J=1,NHEL),K=1,3) /
CLL
     .                   -1,+1,   -1,-1,   -1,+1,
     .                   -1,+1,   -1,-1,   +1,+1,
     .                   -1,+1,   +1,-1,   -1,+1,
     .                   -1,+1,   +1,-1,   +1,+1,
CRR
     .                   -1,+1,   -1,+1,   -1,-1,
     .                   -1,+1,   -1,+1,   +1,-1,
     .                   -1,+1,   +1,+1,   -1,-1,
     .                   -1,+1,   +1,+1,   +1,-1,
CLR
     .                   -1,-1,   -1,-1,   -1,-1,
     .                   -1,-1,   -1,-1,   +1,-1,
     .                   -1,-1,   +1,-1,   -1,-1,
     .                   -1,-1,   +1,-1,   +1,-1/
C--
      DATA        ((( IHEL(I,J,K,2),I=1,NP),J=1,NHEL),K=1,3) /
CLL
     .                   +1,-1,   -1,-1,   -1,+1,
     .                   +1,-1,   -1,-1,   +1,+1,
     .                   +1,-1,   +1,-1,   -1,+1,
     .                   +1,-1,   +1,-1,   +1,+1,
CRR
     .                   +1,-1,   -1,+1,   -1,-1,
     .                   +1,-1,   -1,+1,   +1,-1,
     .                   +1,-1,   +1,+1,   -1,-1,
     .                   +1,-1,   +1,+1,   +1,-1,
CRL
     .                   +1,+1,   -1,+1,   -1,+1,
     .                   +1,+1,   -1,+1,   +1,+1,
     .                   +1,+1,   +1,+1,   -1,+1,
     .                   +1,+1,   +1,+1,   +1,+1/
C--
C>>>
CCDBG       INTEGER*4  JDP(2,NP)
CCDBG       REAL   *4  PVTEST(0:3,10)
C>>>
      DATA NCALL / 0 /
      EXTERNAL   XHFILL
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         AMXM   = SFM(IHANDM+1)
         GMXM   = GMSF(IHANDM+1)
         AMGXM  = AMXM*GMXM
         AMXM2  = AMXM*AMXM
         AMGXM2 = AMGXM*AMGXM
C--
         AMXP   = SFM(IHANDP+1)
         GMXP   = GMSF(IHANDP+1)
         AMGXP  = AMXP*GMXP
         AMXP2  = AMXP*AMXP
         AMGXP2 = AMGXP*AMGXP
C--
         AME   = AMF(1,2,1)
         xALF0 = 1/137.0359895E0
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
C      ( 5) : m(X-)**2
C      ( 6) : m(X+)**2
C      ( 7) : cos(theta_X-)
C      ( 8) : phi_X-
C      ( 9) : cos(theta_l-)    in X- rest frame
C      (10) : phi_l-           in X- rest frame
C      (11) : cos(theta_l+)    in X+ rest frame
C      (12) : phi_l+           in X+ rest frame
C      (13) : e- beam energy spread
C      (14) : e+ beam energy spread
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
cc         CALL PROBEB(REAL(Z(13)),REAL(Z(1)),EBM,EMBM,EMBM0,PEM)
cc         CALL PROBEB(REAL(Z(14)),REAL(Z(2)),EBM,EPBM,EPBM0,PEP)

         call jsfbeamgenr(z(13), z(14), z(1), z(2), ebm, 
     >      embm, epbm, embm0, epbm0, dpdebm)

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
c         DPDEBM    = PEM*PEP
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
         IZE = Z(14)*1.E4
         IF ( MOD(IZE,2).EQ.0 ) THEN
            EMBM = (1-ZG)*EMBM
         ELSE
            EPBM = (1-ZG)*EPBM
         ENDIF
      ENDIF
C--      
      S      = RS*RS
C--
C  Set initial and final states.
C--
      IDPTCL(1, 4) = IGENSF
      IDPTCL(2, 4) = 2
      IDPTCL(3, 4) = 1
      IDPTCL(1, 6) = IGENSF
      IDPTCL(2, 6) = 2
      IDPTCL(3, 6) = 1
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( IHANDM.EQ.IHANDP ) THEN
         IF ( Z(4).LT.HLM ) THEN
            KCMB = 1
            JCMB = IHANDM
            ICMB = NHEL*Z(4)/HLM + 1
         ELSE
            KCMB = 2
            JCMB = IHANDM
            ICMB = NHEL*(Z(4)-HLM)/(1-HLM) + 1
         ENDIF
         WAT  = WAT*NHEL
      ELSE
         KCMB = IHANDM
         JCMB = 3
         ICMB = NHEL*Z(4) + 1
         WAT  = WAT*NHEL
         IF ( KCMB.EQ.1 ) THEN
            WAT = WAT*HLM
         ELSE
            WAT = WAT*(1-HLM)
         ENDIF
      ENDIF
      ICMB = MIN(ICMB,NHEL)
      CALL UVCOPY(NP,IHEL(1,ICMB,JCMB,KCMB),IHLCMB(1))
C--
C  Set invariant masses squared.
C--
      AM3   = SZM(1)
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = SZM(1)
      AM6   = AMF(IDPTCL(1,6),IDPTCL(2,6),IDPTCL(3,6))
C-- X-.
      Q3MN  = AM3 + AM4
      Q3MX  = RS - (AM5+AM6)
      TH3MN = ATAN((Q3MN-AMXM)*(Q3MN+AMXM)/AMGXM)
      TH3MX = ATAN((Q3MX-AMXM)*(Q3MX+AMXM)/AMGXM)
      TH3   = TH3MN + (TH3MX-TH3MN)*Z(5)
      Q32   = AMGXM*TAN(TH3) + AMXM2
      Q3    = SQRT(Q32)
      WAT   = WAT*(TH3MX-TH3MN)*( (Q32-AMXM2)**2 + AMGXM2 )/AMGXM
C-- X+.
      Q4MN  = AM5 + AM6
      Q4MX  = RS - Q3
      TH4MN = ATAN((Q4MN-AMXP)*(Q4MN+AMXP)/AMGXP)
      TH4MX = ATAN((Q4MX-AMXP)*(Q4MX+AMXP)/AMGXP)
      TH4   = TH4MN + (TH4MX-TH4MN)*Z(6)
      Q42   = AMGXP*TAN(TH4) + AMXP2
      Q4    = SQRT(Q42)
      WAT   = WAT*(TH4MX-TH4MN)*( (Q42-AMXP2)**2 + AMGXP2 )/AMGXP
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q32
      XPHASE(1,3) = Q42
C--
C  Angular variables.
C--
      DO 400 IB = 1, 3
         XPHASE(2,IB) = -1 + 2*Z(2*IB+5)
         XPHASE(3,IB) = x2PI  *Z(2*IB+6)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Calculate differential cross sections.
C--
C>>>
CCDBG       PRINT *, '      '
CCDBG       PRINT *, '*************************************'
CCDBG       PRINT *, '    FUNC now calls SGSEEF. '
CCDBG       PRINT *, '*************************************'
CCDBG       PRINT *, '      '
C>>>
      CALL SGSEEF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG      = SG*QED
      FUNC    = SG*WAT
      DSGMDX  = SG
C--      
      IF ( FUNC.LE.0.D0 ) RETURN
C>>>
CCDBG       PRINT *, '      '
CCDBG       PRINT *, ' SG = ', SG,' FUNC = ', FUNC
CCDBG       DO 665 K = 1, NP
CCDBG          JDP(1,K) = 100*IDPTCL(1,K)+10*IDPTCL(2,K)+IDPTCL(3,K)
CCDBG          JDP(2,K) = IHLCMB(K)
CCDBG 665   CONTINUE
CCDBG       PRINT *, ' -------    '
CCDBG       PRINT *, ' IDP, IHL = ', ((JDP(K,KK),K=1,2),KK=1,NP)
CCDBG       DO 666 KK = 1, 3
CCDBG          PRINT *, ' XPHASE(*,', KK, ') = ', (XPHASE(K,KK),K=1,3)
CCDBG 666   CONTINUE
CCDBG       CALL UVZERO(4,PVTEST)
CCDBG       DO 667 KK = 1, NP
CCDBG          PRINT *, ' PV(*,', KK, ') = ', (PV(K,KK),K=0,3)
CCDBG          IF ( KK.GE.3 )
CCDBG      .         CALL PSUMxx(PV(0,KK),PVTEST(0,1),PVTEST(0,1))
CCDBG 667   CONTINUE
CCDBG       PRINT *, ' PVSUM = ', (PVTEST(K,1),K=0,3)
CCDBG       CALL PSUMxx(PV(0,3),PV(0,4),PVTEST(0,3))
CCDBG       CALL PSUMxx(PV(0,5),PV(0,6),PVTEST(0,4))
CCDBG       CALL RDOTxx(PVTEST(0,3),PVTEST(0,3),AM)
CCDBG       AM = SIGN(SQRT(ABS(AM)),AM)
CCDBG       PRINT *, ' AM, PV(S-) = ', AM, (PVTEST(K,3),K=0,3)
CCDBG       CALL RDOTxx(PVTEST(0,4),PVTEST(0,4),AM)
CCDBG       AM = SIGN(SQRT(ABS(AM)),AM)
CCDBG       PRINT *, ' AM, PV(S+) = ', AM, (PVTEST(K,4),K=0,3)
C>>>
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
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
      CALL UVCOPY(4*NP,PVOUT(0,1),PV(0,1))
C--
      ELM   = PV(0,4)
      ELP   = PV(0,6)
      CALL PSUMxx(PV(0,3),PV(0,5),QV(0,1))
      PTOP  = SQRT( QV(1,1)**2 + QV(2,1)**2 )
      COSLM =  PV(3,4)/UABSV3(PV(1,4))
      COSLP = -PV(3,6)/UABSV3(PV(1,6))
      COSLA = (COSLM+COSLP)/2
      PT4   = MAX(SQRT( PV(1,4)**2 + PV(2,4)**2 ),.1)
      PT6   = MAX(SQRT( PV(1,6)**2 + PV(2,6)**2 ),.1)
      ACOP  = ( PV(1,4)*PV(1,6) + PV(2,4)*PV(2,6) )/PT4/PT6
      IF ( ABS(ACOP).GT.1. ) ACOP = SIGN(1.,ACOP)
      ACOP  = xPI - ACOS(ACOP)
      COSSE = XPHASE(2,1)
      PHISE = XPHASE(3,1)
      COSEM = XPHASE(2,2)
      PHIEM = XPHASE(3,2)
      COSEP = XPHASE(2,3)
      PHIEP = XPHASE(3,3)
      CALL XHFILL('h01', DBLE(COSSE)       ,FUNC )
      CALL XHFILL('h02', DBLE(PHISE)*xRD2DG,FUNC )
      CALL XHFILL('h03', DBLE(COSEM)       ,FUNC )
      CALL XHFILL('h04', DBLE(PHIEM)*xRD2DG,FUNC )
      CALL XHFILL('h05', DBLE(COSEP)       ,FUNC )
      CALL XHFILL('h06', DBLE(PHIEP)*xRD2DG,FUNC )
      CALL XHFILL('h07', DBLE(SQRT(Q32))   ,FUNC )
      CALL XHFILL('h08', DBLE(SQRT(Q42))   ,FUNC )
      CALL XHFILL('h09', DBLE(PTOP)        ,FUNC )
      CALL XHFILL('h10', DBLE(ACOP)*xRD2DG ,FUNC )
      CALL XHFILL('h11', DBLE(ELM)         ,FUNC )
      CALL XHFILL('h12', DBLE(ELP)         ,FUNC )
      CALL XHFILL('h13', DBLE(COSLM)       ,FUNC )
      CALL XHFILL('h14', DBLE(COSLP)       ,FUNC )
      CALL XHFILL('h15', DBLE(COSLA)       ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
