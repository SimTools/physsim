C***********************************************************************
C*
C*====================
C* Subroutine PRUJIN
C*====================
C*
C* (Function)
C*    Job initialization for user analysis.
C* (Input & Output )
C*    none.
C* (Relation)
C*    Called from FLANLI
C* (Update Record)
C*    5/29/87  K.Fujii             Original version.
C*
C***********************************************************************
 
      SUBROUTINE PRUJIN
 
      PARAMETER     ( MXxHPK = 400000 )
      COMMON /PAWC/ M(MXxHPK)
C
C========< Entry Point >================================================
C
C--
C  Write user production version #.
C--
      CALL PRPVRS( 'USER:', 001, 900423)
C--
C  Initialize plot buffer.
C--
      CALL HINIT(MXxHPK)
C--
C  Choose appropriate options.
C--
      CALL DOPT('NCALLS',.TRUE.)
      CALL DOPT('L2ZAUT',.TRUE.)
      CALL DOPT('L2COL',2)
C     CALL DOPT('L2COL',1)
C     CALL DOPT('L2SYM',16)
      CALL DOPT('YZERO',.TRUE.)
      CALL HOPTN('STAT',.TRUE.)
      CALL HOPTN('HOMIT',.TRUE.)
C--
C  Initialize hists and plots.
C--
      CALL HDEF1(  1,'I*4',220,  0.,0.005   ,  'x-  @')
      CALL HDEF1(  2,'I*4',220,  0.,0.005   ,  'x+  @')
      CALL HDEF1(  3,'I*4',220,  0.,0.005   ,  'rsh @')
C--
C  2-d plots.
C--
      CALL HDEF2(201,'I*4',210,210,0.,0.,.005,.005,'(x-,x+) @')
C--
C  That's it.
C--
      RETURN
      END
