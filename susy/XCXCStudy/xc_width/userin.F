C* (Update Record)
C*   95/04/10  K.Fujii		Original version for X+ decay.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
      PARAMETER (MXDIM = 50 )
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2
C--
#include "XCXCCN.inc"
#include "USRPRM.inc"
#include "USSPRM.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*4        ALF
C
C========< Entry Point >================================================
C
C--
C  Check if inputs from FT05 are consistent with those used in the
C  integration step, when NGETBS = 1. Otherwise, Write the inputs 
C  to 'BASES:Inputs'.
C--
      CALL SPPARM
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
      CALL INSSCN(ALF,ALFS,AMSW,AMSZ,AM0,AMU,AM2,TANB,AMA,ASFT)
C--
C  Set BASES constants.
C--

      NDIM   = 7
      NWILD  = 4
      IOFF   = 0
C>>>
C      NOIG = NWILD + 1
      NOIG = 999
C>>>
C--
      NZZ = NDIM
C--
      NCALL  = NCAL
      ITMX1  = ITM1
      ITMX2  = ITM2    
C--
C  Set lower and upper limits.
C--
      DO 100 IX = 1, NDIM
         XL(IX) = 0
         XU(IX) = 1
         IF ( IX.GE.NOIG ) IG(IX) = 0
         ISHUFL(IX) = IX + IOFF
100   CONTINUE
C--
C  Define hists. and plots.
C--
      QMX  = 250
      CALL XHINIT( 1, -1.0D0, 1.0D0, 50,'cos(theta_a+)      ')
      CALL XHINIT( 2,  0.0D0,360.D0, 50,'phi_b+             ')
      CALL XHINIT( 3,  0.0D0,   QMX, 50,'M_ab+              ')
      CALL XHINIT( 4,  0.0D0,   QMX, 50,'M_ac+              ')
      CALL XHINIT( 5,  1.0D0,  9.D0,  8,'Hel. comb.         ')
      CALL XHINIT( 6,  0.0D0,   QMX, 50,'E_e/mu             ')
      CALL XHINIT( 7, -1.0D0, 1.0D0, 50,'cos(theta_l-)      ')
      CALL XHINIT( 8,  1.0D0, 13.D0, 12,'Decay mode( X+ )   ')
C--
C  That's it.
C--
      RETURN
      END
