CC**********************************************************************
C*
C*===================
C* Subroutine FLANLI
C*===================
C*
C*(Purpose)
C*    Define analysis program. This version is for BASES5.1.
C*(Relation)
C*    Called by PRMAIN.
C*(Author)
C*   17-May-1990  A.Miyamoto	For Full system
C*   15-Mar-1995  K.Fujii	Added DSDBGB, DSRCOV, and DSSQUZ.
C*
CC**********************************************************************
C
      SUBROUTINE FLANLI
C
#include "fldscr.inc"
C
C======< ENTRY POINT >==================================================
C
C ----------------------------------------------------------------------
C (1) Declare version number of Full module.
C ----------------------------------------------------------------------
C
      CALL PRPDSC( 'JLC-BEAM', 101, 950623)
C
C ----------------------------------------------------------------------
C (2) Reset counters.
C ----------------------------------------------------------------------
C
      NPROG = 0
      NENTR = 0
      NHIST = 0
C
C ----------------------------------------------------------------------
C (3) Define analysis routine and calling sequence.
C ----------------------------------------------------------------------
C
C                     Define BASES generator.
      CALL DSBS51

C                     Define User routine.
      CALL DSCUSR
 
      RETURN
      END

