C*
C* (Update Record)
C*   95/Oct/14  T.Tsukamoto		Modified for BASES25.
C*   96/Mar/08  K.Fujii			Modified for BASES5.1.
C*   99/05/20   K.Fujii			Modified for JFS.
C*
      SUBROUTINE USERIN
 
C ....................................................................
C
C          E(+) E(-) INTO E(+) E(-) MU(+) MU(-)
C                                                        20 JAN. 1986
C
C          EXACT CALCULATION
C          CROSS SECTION IN NB =10**(-33) CM**2
C
C......................................................................
 
      IMPLICIT REAL*8(A-H,O-Z)
C--
C  User's commons.
C--
#include "USMPRM.inc"
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "USRPRM.inc"
#include "KINEM.inc"
#include "CUTS.inc"
C--
C  Local arrays.
C--
      REAL   *8    QF(3,2,2)
C
C========< Entry Point >================================================
C
C--
C  Set fermion charge: QF(IG,IT,LQ), where
C     IG = generation number
C     IT = (1,2) = (up-type,down-type)
C     LQ = (1,2) = (lepton,quark)
C--
      DO 10 IGN = 1, 3
         QF(IGN,1,1) =  0
         QF(IGN,2,1) = -1
         QF(IGN,1,2) =  2.D0/3
         QF(IGN,2,2) = -1.D0/3
10    CONTINUE
C--
C  Initialize constants. 
C--
      ALF    = 1/ALFI
      CALL INSMCN(ALF,DBLE(ALFS),DBLE(AMSW),DBLE(AMSZ),
     .                           DBLE(AMSH),DBLE(AMST))
C--
      CE     = SQRT(4*xPI*ALF)
      EM     = AMF(1,2,1)
C--
C  Set the charge and mass of the final-state fermions.
C--
      CHARGE = QF (IGFR,ITFR,LQFR)
      UM     = AMF(IGFR,ITFR,LQFR)
C--
C  Center of mass energy.
C--
      W      = SQRTS
      E      = W/2
      S      = W*W
C--
C  Derived quantities.
C--
      CE8    = CE**8
C--
      EME    = EM/E
      EME2   = EME**2
      EME4   = EME2**2
      EME6   = EME2*EME4
      EME8   = EME**8
C--
      UME    = UM/E
      UME2   = UME**2
      UME4   = UME2**2
      UME6   = UME2*UME4
      UME8   = UME4**2
C--
      BETHA  = SQRT((1+EME)*(1-EME))
      FACT   = CE8/(4*E**2*(2*xPI)**7)*xGV2PB/BETHA
C--
C  Cuts.
C--
C e+
      PCUT    = PECT
      PCUT    = PCUT/E
      CTAGL_2 = CECL
      CTAGU_2 = CECU
C e-
      PTAGC   = PPCT
      PTAGC   = PTAGC/E
      CTAGL_1 = CPCL
      CTAGU_1 = CPCU
C f
      PMUCUT  = PFCT
      PMUCUT  = PMUCUT/E
      CMUCUT  = CFCT
      WMIN    = WMNF
      WMIN    = WMIN/E
C--
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     W      = ', W, ' GeV'
C      PRINT *, '     SGEBM  = ', SGEBM
C      PRINT *, '     POLE   = ', POLE
C      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     CHARGE = ', CHARGE
      PRINT *, '     UM     = ', UM, ' GeV'
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
       RETURN
       END




