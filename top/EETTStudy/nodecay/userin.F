C* (Update Record)
C*   95/04/29  K.Fujii		Original version for EETT.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
      PARAMETER (MXDIM = 50 )
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2
C--
#include "EETTCNF.inc"
#include "USRPRM.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
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
      CALL INSMCN(ALF,DBLE(ALFS),DBLE(AMSW),DBLE(AMSZ),
     .                           DBLE(AMSH),DBLE(AMST))
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
      IMD1MN = IMD1LO
      IMD1MX = IMD1HI
      IMD2MN = IMD2LO
      IMD2MX = IMD2HI
C--
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     IMD1MN = ', IMD1MN
      PRINT *, '     IMD1MX = ', IMD1MX      
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX      
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  Set BASES constants.
C--
      IF ( ISRBM.EQ.1 ) THEN
         NDIM   = 9
         NWILD  = 5
         IOFF   = 3
      ELSE IF ( ISRBM.EQ.2 ) THEN
         NDIM   = 10
         NWILD  = 6
         IOFF   = 2
      ELSE IF ( ISRBM.EQ.3 ) THEN
         NDIM   = 14
         NWILD  = 8
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
      CALL XHINIT( 1,   0.D0, 1.0D0, 50,'E_1/E_bm      ')
      CALL XHINIT( 2, -1.0D0, 1.0D0, 50,'cos_1         ')
      CALL XHINIT( 3,  0.0D0,360.D0, 50,'phi_1         ')
      CALL XHINIT( 4,   0.D0, 1.0D0, 50,'E_2/E_bm      ')
      CALL XHINIT( 5, -1.0D0, 1.0D0, 50,'cos_2         ')
      CALL XHINIT( 6,  0.0D0,360.D0, 50,'phi_2         ')
      CALL XHINIT( 7,   0.D0, 1.0D0, 50,'TT/ROOTS      ')
      CALL XHINIT( 8, -1.0D0, 1.0D0, 50,'cos_t         ')
      CALL XHINIT( 9,  0.0D0,360.D0, 50,'phi_t         ')
      CALL XHINIT(10,  0.0D0, 1.0D0, 50,'RS/ROOTS      ')
      CALL XHINIT(11,  1.0D0, 9.0D0,  8,'helicity      ')
      CALL XHINIT(12, -20.D0,20.0D0, 50,'eta_1         ')
      CALL XHINIT(13, -20.D0,20.0D0, 50,'eta_2         ')
      CALL XHINIT(14,   0.D0, 1.0D0, 50,'PT_t/E_bm     ')
      CALL XHINIT(15,   0.D0, 1.0D0, 50,'PT_tb/E_bm    ')
      CALL XHINIT(16,   0.D0, 1.0D0, 50,'E_1/E_bm      ')
      CALL XHINIT(17,   0.D0, 1.0D0, 50,'E_2/E_bm      ')
      CALL XHINIT(18,   0.D0, 1.0D0, 50,'E_t/E_bm      ')
      CALL XHINIT(19,   0.D0, 1.0D0, 50,'E_tb/E_bm     ')
C--
C  That's it.
C--
      RETURN
      END
