C***********************************************************************
C*
C*======================================----==
C* Subroutine PRUSER(IDREC,LEVEL,IDEBUG,IRET)
C*======================================----==
C*
C* (Function)
C*     User Event analysis routine.
C* (Input)
C*     IDREC   :  Record ID.
C*     ILEVEL  :  Level of analysis.
C*     IDEBUG  :  Debug flag.
C* (Output)
C*     IRET    :  Return flag.
C* (Update Record)
C*   95/06/15  K.Fujii             Original version.
C*
C***********************************************************************
 
      SUBROUTINE PRUSER( IDREC, LEVEL, IDEBUG, IRET )
 
      INCLUDE 'BEMCNF.inc'
C
C==================<<  Entry Point  >>================================
C
C--
C  Reset return code, Skip except EVENT record.
C--
      IRET = 0
      IF( IDREC .NE. 1 ) RETURN
      LEVEL  = LEVEL
      IDEBUG = IDEBUG
      IRET = -1
C--
C  Hists. and plots here.
C--
      CALL HCUM1(  1,XM,1.)
      CALL HCUM1(  2,XP,1.)
      CALL HCUM1(  3,FS,1.)
      CALL HCUM2(201,XM,XP,1.) 
C--
C  That's it.
C--
      IRET = 0
9999  RETURN
      END
