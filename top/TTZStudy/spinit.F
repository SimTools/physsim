CC**********************************************************************
C*===================
C* Subroutine SPINIT
C*===================
C*
C*     This routine is called at Run initialization.
C*     Store Run conditions into Begin_Run records.
C*
C* (Update Record)
C*   90/05/17  A.Miyamoto	Original version.
C*   95/03/15  K.Fujii		Deleted SPFUNC call for BASES V5.1.
C*
C*
CC**********************************************************************
C
      SUBROUTINE SPINIT
C
#include "TTZCNF.inc"
      COMMON /SSWORK/ ELMNO2(20), ELMNO3(10)
      CHARACTER*32    PRCNAM
C
C
C =====< Entry Point >=================================================
C
C ---------------------------------------------------------------------
C (1) Store Run parapemter into Bank 'Spring' Elm# 2.
C ---------------------------------------------------------------------
C
C     ELMNO2(1) = E-beam (GeV)
C           (2) = 1. When Beamstrahlung is on
C               = 0. When Beamstrahlung is off
C           (3:12) = reserved to store beam parameter.
C
      ELMNO2(1) = ROOTS/2
      ELMNO2(2) = 0.
      ELMNO2(3) = 0.
      ELMNO2(4) = 0.
      ELMNO2(5) = 0.
      ELMNO2(6) = 0.
      ELMNO2(7) = 0.
      ELMNO2(8) = 0
      ELMNO2(9) = 0
      ELMNO2(10)= 0
      ELMNO2(11)= 0
      ELMNO2(12)= 0
C
      PRCNAM    = 't-tbar-Z'
      READ(PRCNAM,'(8A4)') (ELMNO2(K),K=13,20)
      NWD       = 20
C
      CALL TBPUT(2,'Spring',2,NWD,ELMNO2, IRET )
      IF ( IRET .LT. 0 ) THEN
        PRINT *,' USERIN failed to TBPUT SPRING Elm#2 ',IRET
        STOP
      ENDIF
C
C ---------------------------------------------------------------------
C (2) Store Run parapemter into Bank 'Spring' Elm# 3.
C     Store information specific to this process.
C ---------------------------------------------------------------------
C
      RETURN
      END
