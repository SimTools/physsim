CC**********************************************************************
C*
C*====================
C* Subroutine PRUCOM
C*====================
C*
C* (Function)
C*   Read parameters for PRUSER.
C*   This routine is called prior to the job initialization by the
C*   steering routine when the string 'USER' is read from FT05.
C*   This routine should return to its caller when read 'EXIT'.
C* (Update Record)
C*   95/03/17  K.Fujii		Original version for analysis.
C*
CC**********************************************************************

      SUBROUTINE PRUCOM

      INCLUDE 'OUTHBK.inc'
      CHARACTER*4   CMD
C
C =====< Entry Point >==================================================
C
C --------------------------------------------
C  Example of read parameters using COM/CLI
C --------------------------------------------
C
      PRINT *, 'USER Analysis'
      PRINT *,' '
C--
C  Default job parameters.
C--
      HFILEN = 'ttanl.hbook'
C--
C  Loop over input cards.
C--
100   CONTINUE
      PRINT *, 'Parameter keys for US: DISF'
C--
C  Read a command key.
C--
      CMD = '?'
      CALL TXTRD ( ' Enter parameter name or EXIT ', CMD )
      CALL STRUPC( 4, CMD )
C--
C  Branch on keys.
C--
      IF ( CMD(1:4).EQ.'EXIT') THEN
	 PRINT *,'End of parameter setting (in PRUCOM).'
	 					GO TO 900
C--
C  USER I/O parameters.
C--
      ELSE IF ( CMD(1:4).EQ.'DISF') THEN
	 HFILEN = '?'
	 CALL TXTRD('DISF: ',HFILEN)
C--
C  USER cuts, etc.
C--
C      ELSE IF ( CMD(1:4).EQ.'ACUT') THEN
C	 CALL FLTRD('ACUT: ',ACUT)
C      ELSE
C--
C  Help for USER inputs.
C--
	 PRINT *,' '
	 PRINT *,'Command parameters for USER:'
	 PRINT *,'DISF: output disfile name.'
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
C  Parameters for PRUSER.
C--
      PRINT *,' '
      PRINT *,'Parameters for USER'
      PRINT *,' '
      PRINT *,' DISF=', HFILEN(:LNBLNK(HFILEN))
C--
C  That's it.
C--
      RETURN
      END
