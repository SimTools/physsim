C* (Update Record)
C*  95/06/23  K.Fujii	Original version for BASES 5.1.
C*
      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
      PARAMETER (MXDIM = 50 )
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2
C--
      INCLUDE 'BEMCNF.inc'
      INCLUDE 'USRPRM.inc'
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
C  Set user constants.
C--
      ROOTS  = SQRTS
      SGEBM  = SGMEBM
      ISRBM  = ISRB
      IBMTP  = IBTP
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     IBMTP  = ', IBMTP
      PRINT *, '  '
      PRINT *, ' ********************************** '
      IF ( ISRBM.EQ.1 ) THEN
         NDIM   = 1
         NWILD  = 1
         IOFF   = 4
      ELSE IF ( ISRBM.EQ.2 ) THEN
         NDIM   = 4
         NWILD  = 4
         IOFF   = 0
      ELSE IF ( ISRBM.EQ.3 ) THEN
         NDIM   = 5
         NWILD  = 5
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
C  Set lower and upper limits.
C--
      CALL XHINIT( 1,0.0D0, 3.0D0, 50,'sqrt(s)_1    ')
      CALL XHINIT( 2,3.0D0, 4.5D0, 50,'sqrt(s)_2    ')
      CALL XHINIT( 3,4.5D0, 6.0D0, 50,'sqrt(s)_3    ')
      CALL XHINIT( 4,6.0D0, 7.5D0, 50,'sqrt(s)_3    ')
      CALL XHINIT( 5,0.5D0, 1.0D0, 50,'rsh/roots    ')
C--
C  That's it.
C--
      RETURN
      END
