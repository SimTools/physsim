CC**********************************************************************
C*
C*=========================================================
C* Subroutine AMPNN(GWF, AMW,GMW, SPEM,SPEP, SP1,SP2, AMP)
C*====================================================---==
C*
C* (Purpose)
C*    Calculate amplitude for e- e+ ---> nu nubar.
C*    Inputs are spinors thus replaceable with spinors made from 
C*    final stable particles.
C* (Inputs)
C*       GWF      : (R*8) : weak coupling.
C*       AMW      : (R*8) : W mass.
C*       GMW      : (R*8) : W width.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state nu spinor.
C*       SP2 (*)  : (C*16) : final-state nub spinor.
C* (Output)
C*       AMP      : (C*16) : total amplitudes.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    2005/03/04  K.Fujii	Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPNN(GWF, AMW,GMW, SPEM,SPEP, SP1,SP2, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GWF(2), AMW, GMW
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), AMP
C--
      COMPLEX*16 WRK(6)
C--
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Calculate internal lines.
C--
      CALL JIOXXX(SPEM,SP1 ,      GWF, AMW,GMW, WRK)
      CALL IOVXXX(SP2 ,SPEP, WRK, GWF, AMP)
C--
C  That's it.
C--
      RETURN
      END
