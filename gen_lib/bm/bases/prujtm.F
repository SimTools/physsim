C***********************************************************************
C*
C*====================
C* Subroutine PRUJTM
C*====================
C*
C* (Purpose)
C*    Job termination for user analysis.
C* (Input & Output )
C*    none.
C* (Relation)
C*    Called by FLANLI
C* (Update Record)
C*   95/03/06  K.Fujii		Orignal version.
C*
C***********************************************************************
 
      SUBROUTINE PRUJTM

      INCLUDE 'OUTHBK.inc'
C
C========< Entry Point >================================================
C
C--
C  Output hists. and plots.
C--
      CALL HRPUT(0,HFILEN,'N')
C      CALL HPRINT(0)
C--
C  That's it.
C--
      RETURN
      END
