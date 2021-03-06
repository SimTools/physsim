CC**********************************************************************
C*
C*====================
C* Subroutine BS51CM
C*====================
C*
C* (Function)
C*   Read parameters for a BASES/SPRING job.
C*   This routine is called prior to the job initialization by the
C*   steering routine when the string 'BASES' is read from FT05.
C*   This routine should return to its caller when read 'EXIT'.
C* (Update Record)
C*   94/10/19  A.Miyamoto       Derived from SPGPAR.
C*   95/06/23  K.Fujii		For BM generator.
C*
CC**********************************************************************

      SUBROUTINE BS51CM

      INCLUDE 'USRPRM.inc'
#include "BS51PR.inc"
#include "BSITPR.inc"
      CHARACTER*4   CMD
C
C =====< Entry Point >==================================================
C
C --------------------------------------------
C  Example of read parameters using COM/CLI
C --------------------------------------------
C
      PRINT *, 'BASES_V5.1 - Change parameters for BASES V5.1 ',
     .         ' generator .'
      PRINT *,' '
C--
C  Default user parameters.
C--
      BSHSTF = 'bases.tdr'
      LUHIST = 24
C--
      NCAL   = 400000
      ITM1   = 5
      ITM2   = 10
C--
C  Set default beam parameters.
C--
      SQRTS  = 300.
      SGMEBM = 0.005
      ISRB   = 2
      IBTP   = 33
C--
C  Loop over input cards.
C--
100   CONTINUE
      PRINT *, 'Parameter keys for BS: GETB PUTB DOSP PHST PINF PLOG ',
     .         'INF OUTF HSTF'
      PRINT *, '                     : NCAL ITM1 ITM2'
      PRINT *, '               for SM: ALFI ALFS AMSW AMSZ AMSH AMST '
      PRINT *, '               for BM: ECM  SGEB ISRB IBTP'
C--
C  Read a command key.
C--
      CMD = '?'
      CALL TXTRD ( 'Enter parameter name or EXIT ', CMD )
      CALL STRUPC( 4, CMD )
C--
C  Branch on keys.
C--
      IF ( CMD(1:4).EQ.'EXIT') THEN
	 PRINT *,'End of parameter setting (in BS51CM).'
	 					GO TO 900
C--
C  BASES V5.1 parameters.
C--
      ELSE IF ( CMD(1:4).EQ.'GETB') THEN
	 CALL INTRD('GETB: ',NGETBS)
      ELSE IF ( CMD(1:4).EQ.'PUTB') THEN
	 CALL INTRD('PUTB: ',NPUTBS)
      ELSE IF ( CMD(1:4).EQ.'DOSP') THEN
	 CALL INTRD('DOSP: ',NDOSPR)
      ELSE IF ( CMD(1:4).EQ.'PHST' ) THEN
	 CALL INTRD('PHST: ',NPRHST)
      ELSE IF ( CMD(1:4).EQ.'PINF') THEN
	 CALL INTRD('PINF: ', NPRINF)
      ELSE IF ( CMD(1:4).EQ.'PLOG') THEN
	 CALL INTRD('PLOG: ',NBSINF)
      ELSE IF ( CMD(1:3).EQ.'INF') THEN
	 BSINF = ' '
	 CALL TXTRD('Input File: ',BSINF)
      ELSE IF ( CMD(1:4).EQ.'OUTF') THEN
	 BSOUTF = ' '
	 CALL TXTRD('Output File: ',BSOUTF)
      ELSE IF ( CMD(1:4).EQ.'HSTF') THEN
	 BSHSTF = ' '
	 CALL TXTRD('Histgram File: ',BSHSTF)
C--
      ELSE IF ( CMD(1:4).EQ.'NCAL') THEN
	 CALL INTRD('NCAL: ',NCAL)
      ELSE IF ( CMD(1:4).EQ.'ITM1') THEN
	 CALL INTRD('ITM1: ',ITM1)
      ELSE IF ( CMD(1:4).EQ.'ITM2') THEN
	 CALL INTRD('ITM2: ',ITM2)
C--
C  BM parameters.
C--
      ELSE IF ( CMD(1:3).EQ.'ECM') THEN
	 CALL FLTRD('ECM: ',SQRTS)
      ELSE IF ( CMD(1:4).EQ.'SGEB') THEN
	 CALL FLTRD('SGEB: ',SGMEBM)
      ELSE IF ( CMD(1:4).EQ.'ISRB') THEN
	 CALL INTRD('ISRB: ',ISRB)
      ELSE IF ( CMD(1:4).EQ.'IBTP') THEN
	 CALL INTRD('ISRB: ',IBTP)
      ELSE
C--
C  Help for BASES inputs.
C--
	 PRINT *,' '
	 PRINT *,'Command parameters for BASES_V5.1:'
	 PRINT *,'GETB: 0 = Do BASES calculation at begin run.'
	 PRINT *,'      1 = Get BASES data from TBS bank.'
         PRINT *,'      2 = Get BASES data from BASES ',
     .           'standard format file.' 
	 PRINT *,'PUTB: 0 = Does not output BASES result.'
	 PRINT *,'      1 = Write BASES results into TBS bank.'
	 PRINT *,'      2 = Write BASES results with BASES ',
     .           ' standard file format.'
	 PRINT *,'DOSP: 1(0) = Do (not do) SPRING.'
	 PRINT *,'PHST: 1(0) = Print (not print) BASES/SPRING  Hist.'
	 PRINT *,'PINF: 1(0) = Print (not print) BASES/SPRING Info.'
	 PRINT *,'PLOG: 1(0) = Print (not print) BASES results.'
	 PRINT *,'INF:  Input file name for GETB=2 case.'
	 PRINT *,'OUTF: Output file name for PUTB=2 case.'
	 PRINT *,'HSTF: Histogram file name.'
	 PRINT *,'Valid combination of parameters are:'
	 PRINT *,' (1)  GETB=0, PUTB=0/1/2, DOSP=0/1 '
         PRINT *,' (2)  GETB=1, PUTB=0/2, DOSP=1 '
	 PRINT *,' (3)  GETB=2, PUTB=0/1/2, DOSP=1/2'
      ENDIF
C--
C  Go to next card.
C--
      						GO TO 100
C--
C     Call FLTRD( 'message', XDATA ) ; Read REAL*4 value.
C     Call INTRD( 'message', IDATA ) ; Read Integer*4 calue.
C     Call TXTRD( 'message', STRING ) ; Read String
C     Call TITRD( 'message', STRING ) ; Read a line.
C--
C
C --------------------------------------------
C  Print out the modified parameters.
C --------------------------------------------
C
900   CONTINUE
C--
C  Parameters for BASES V5.1.
C--
      PRINT *,' '
      PRINT *,'Parameters for BASES V5.1'
      PRINT *,' '
      PRINT *,' GETB=',NGETBS,' PUTB=',NPUTBS,' DOSP=',NDOSPR
      PRINT *,' PHST=',NPRHST,' PINF=',NPRINF,' PLOG=',NBSINF
      IF( NGETBS.EQ.2 ) PRINT *,' In file=',BSINF(:LNBLNK(BSINF))
      IF( NPUTBS.EQ.2 ) PRINT *,' Out file=',BSOUTF(:LNBLNK(BSOUTF))  
      PRINT *,' HSTF=',BSHSTF(:LNBLNK(BSHSTF))
C--
C  User parameters.
C--
      PRINT *,' '
      PRINT *,'User Parameters'
      PRINT *,' '
      PRINT *,' ECM =',SQRTS
      PRINT *,' SGEB=',SGMEBM
      PRINT *,' ISRB=',ISRB
      PRINT *,' IBTP=',IBTP
      PRINT *,' '
C--
C  That's it.
C--
      RETURN
      END
