C*
C* (Update Record)
C* 95/03/15  K.Fujii	Modified for BASES V5.1 in FULL.
C* 95/03/21  K.Fujii	Now calls SPPARM.
C* 95/06/07  K.Fujii	Added ISRB switch.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL FUNC
C--
      PARAMETER (MXDIM = 50 )
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2
C--
#include "USRPRM.inc"
#include "TTZCNF.inc"
#include "BEMCNS.inc"
#include "USMPRM.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*4   ALF
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
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      POLE   = POLEBM
      SGEBM  = SGMEBM
      ISRBM  = ISRB
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
      IF ( ROOTS.LT.AMST+AMST+AMSZ ) THEN
         PRINT *, ' ROOTS is too small. '
         PRINT *, '    AMT = ', AMST
         PRINT *, '    AMZ = ', AMSZ
      ENDIF
C--
C  Set BASES constants.
C--
      IF ( ISRBM.EQ.1 ) THEN
         NDIM   = 7
         NWILD  = 3
         IOFF   = 3
      ELSE IF ( ISRBM.EQ.2 ) THEN
         NDIM   = 8
         NWILD  = 4
         IOFF   = 2
      ELSE IF ( ISRBM.EQ.3 ) THEN
         NDIM   = 12
         NWILD  = 6
         IOFF   = 0
      ELSE
         PRINT *, ' >>> USERIN : Invalid ISRBM = ', ISRBM
         PRINT *, '   Will STOP immediately. '
         STOP
      ENDIF
C>>>
      NOIG = NWILD + 1
C      NOIG = 999
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
      CALL XHINIT( 1, -1.0D0, 1.0D0, 50,'cos(theta_Z) ')
      CALL XHINIT( 2,  0.0D0,360.D0, 50,'phi_Z        ')
C      CALL XHINIT( 3,   0.D0,  1.D0, 50,'m(t-tb)/roots')
      CALL XHINIT( 3, 200.D0,700.D0, 50,'m(t-tb)      ')
      CALL XHINIT( 4, -1.0D0, 1.0D0, 50,'cos(theta_t) ')
      CALL XHINIT( 5,  0.0D0,360.D0, 50,'phi_t        ')
      CALL XHINIT( 6, -1.0D0, 1.0D0, 50,'cos(theta_t)_lab    ')
      CALL XHINIT( 7,  1.0D0, 25.D0, 24,'helicity combination')
      CALL DHINIT( 9,0.D0,1.D0,50,-1.D0,1.D0,50,'E_Z/E_bm-cos(th_Z)')
      CALL DHINIT(10,0.D0,1.D0,50, 0.D0,1.D0,50,'E_t/E_bm-E_tb/E_bm')
C--
C  That's it.
C--
      RETURN
      END
