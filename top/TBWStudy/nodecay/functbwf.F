C* 2010/11/29  K.Fujii	Derived from TTStudy/functtf.F
C*
      REAL*8 FUNCTION FUNC(ZZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
      INCLUDE 'TBWCNF.inc'
#include "SMPTAB.inc"
#include "SMCONS.inc"
#include "BSHUFL.inc"
C--
C  # external lines (NP) and helicity combinations (NHEL).
C--
      PARAMETER  ( NP = 5, NHEL = 6 )
      REAL   *8   SG, PV(0:3,NP), QV(0:3,10)
C--
      INTEGER*4  IHEL(NP,2*NHEL)
      DATA ( ( IHEL(I,J),I=1,NP ), J= 1,2*NHEL )  /
     .             -1,+1,   +1, -1,-1,
     .             -1,+1,   +1, +1,-1,
     .             -1,+1,   +1, -1, 0,
     .             -1,+1,   +1, +1, 0,
     .             -1,+1,   +1, -1,+1,
     .             -1,+1,   +1, +1,+1,
     .             +1,-1,   +1, -1,-1,
     .             +1,-1,   +1, +1,-1,
     .             +1,-1,   +1, -1, 0,
     .             +1,-1,   +1, +1, 0,
     .             +1,-1,   +1, -1,+1,
     .             +1,-1,   +1, +1,+1/

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
         AME   = AMF(1,2,1)
         xALF0 = 1/137.0359895D0
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
C      ( 1) : e- beam
C      ( 2) : e+ beam
C      ( 3) : bremsstrahlng
C      ( 4) : e- helicity
C      ( 5) : final state helicity combination
C      ( 6) : m(t_bar)**2
C      ( 7) : cos(theta_t)
C      ( 8) : phi_t
C      ( 9) : cos(theta_q_bar) in t_bar rest frame
C      (10) : phi_q_bar        in t_bar rest frame
C      (11) : e- beam gaussian spread
C      (12) : e+ beam gaussian spread
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
      IF ( ISRBM.EQ.3 ) THEN
         call jsfbeamgenr(Z(1), Z(2), Z(11), Z(12), ebm,
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
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      IF ( ISRBM.GE.2 ) THEN
         BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
         QED = ( 1 + 3*BTE/4 )*( 1 + 2*xALF0/xPI*( xPI**2/6 - 1/4.D0 ) )
         ZG  = Z(3)**(1/BTE)
         RS  = SQRT(1-ZG)*RS
         IZE = Z(8)*1.D10
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
C     (3,4)=(bbar,t)
C--
      IDPTCL(1, 3) = 3
      IDPTCL(2, 3) = 2
      IDPTCL(3, 3) = 2
      IDPTCL(1, 4) = 3
      IDPTCL(2, 4) = 1
      IDPTCL(3, 4) = 2
C--
C  Select helicity combination.
C--
      HLM   = (1-POLE)/2
      IF ( Z(4).LT.HLM ) THEN
         ICMB = 1
      ELSE
         ICMB = 2
      ENDIF
      JCMB = NHEL*Z(5)+1
      JCMB = MIN(NHEL,JCMB)
      JCMB = (ICMB-1)*NHEL + JCMB
      CALL UVCOPY(NP,IHEL(1,JCMB),IHLCMB(1))
      WAT = WAT*NHEL
C--
C  Set invariant masses squared.
C--
C--
      AM3   = AMF(IDPTCL(1,3),IDPTCL(2,3),IDPTCL(3,3))
      AM4   = AMF(IDPTCL(1,4),IDPTCL(2,4),IDPTCL(3,4))
      AM5   = AMW
C-- t_bar.
      Q3MN  = AM3 + AM5
      Q3MX  = RS - AM4
      Q32MN = Q3MN*Q3MN
      Q32MX = Q3MX*Q3MX
      Q32   = Q32MN + (Q32MX-Q32MN)*Z(6)
      Q3    = SQRT(Q32)
      WAT   = WAT*(Q32MX-Q32MN)
C-- t.
      Q42   = AM4*AM4
      Q4    = AM4
C-- W-.
      Q52   = AM5*AM5
      Q5    = AM5
C--
      XPHASE(1,1) = RS*RS
      XPHASE(1,2) = Q32
C--
C  Angular variables.
C--
      DO 400 IB = 1, 2
         XPHASE(2,IB) = -1 + 2*Z(2*IB+5)
         XPHASE(3,IB) = x2PI  *Z(2*IB+6)
         WAT = WAT*2*x2PI
400   CONTINUE
C--
C  Calculate differential cross sections.
C--
      CALL SGTBWF(IDPTCL,IHLCMB,XPHASE,SG,PV)
      SG     = SG*QED
      FUNC   = WAT*SG
      DSGMDX = SG
C--
C  Boost everything to laboratroy frame.
C--
      QV(0,1) = EMBM + EPBM
      QV(1,1) = 0
      QV(2,1) = 0
      QV(3,1) = EMBM - EPBM
      IF ( ABS(QV(3,1)).GT.1.D-3 ) THEN
C>>>
CCDBG       PRINT *, ' QV = ', (QV(K,1),K=0,3)
C>>>
         DO 500 IP = 1, NP
            CALL BOOSTx(PV(0,IP),QV(0,1),PVOUT(0,IP))
500      CONTINUE
      ELSE
         CALL UVCOPY(2*4*NP,PV(0,1),PVOUT(0,1))
      ENDIF
C--
C  Fill hists. and plots.
C--
      COST  = XPHASE(2,1)
      PHIT  = XPHASE(3,1)
      COSBB = XPHASE(2,2)
      PHIBB = XPHASE(3,2)
      PTOP  = (RS/2)*BETA(Q42/RS**2,Q32/RS**2)
      CALL PSUMxx(PVOUT(0,4),PVOUT(0, 3),QV(0,9))
      CALL PSUMxx(PVOUT(0,3),PVOUT(0, 5),QV(0,3))
      CALL UVCOPY(2*4,PVOUT(0,4),QV(0,4))
      CALL RDOTxx(QV(0,9),QV(0,9),QTBB2)
      QTBB = SQRT(QTBB2)
      PW   = SQRT(PVOUT(1,5)**2+PVOUT(2,5)**2+PVOUT(3,5)**2)
      CSWM = PVOUT(3,5)/PW

      PT1  = SQRT(QV(1,3)**2+QV(2,3)**2+QV(3,3)**2)
      PT2  = SQRT(QV(1,4)**2+QV(2,4)**2+QV(3,4)**2)
C--
      CALL XHFILL( 'h01', COST              ,FUNC )
      CALL XHFILL( 'h02', PHIT*xRD2DG       ,FUNC )
      CALL XHFILL( 'h03', COSBB             ,FUNC )
      CALL XHFILL( 'h04', PHIBB*xRD2DG      ,FUNC )
      CALL XHFILL( 'h07', RS/ROOTS          ,FUNC )
      CALL XHFILL( 'h08', QTBB              ,FUNC )
      CALL XHFILL( 'h09', Q4                ,FUNC )
      CALL XHFILL( 'h10', Q3                ,FUNC )
      CALL XHFILL( 'h11', PTOP              ,FUNC )
      CALL XHFILL( 'h12', PT1               ,FUNC )
      CALL XHFILL( 'h13', PT2               ,FUNC )
      CALL XHFILL( 'h14', JCMB + .1D0       ,FUNC )
      CALL XHFILL( 'h15', CSWM              ,FUNC )
C--
C  That's it.
C--
      RETURN
      END
