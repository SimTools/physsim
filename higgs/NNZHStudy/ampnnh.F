CC**********************************************************************
C*
C*=====================================================
C* Subroutine AMPNNH(GWF, GWWH, AMW,GMW, AME, AMH,GMH,
C*                       SPEM,SPEP, SP1,SP2, SCH, AMP)
C*================================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ e- ---> nu nub h.
C*    Inputs are spinors and a scalar wave functionthus 
C*    replaceable with virtual spinors and a scalar wave function
C*    made from final stable particles.
C* (Inputs)
C*       GWF      : (R*8)  : weak coupling.
C*       GWWH     : (R*8)  : WWH coupling.
C*       AMW      : (R*8)  : W mass.
C*       GMW      : (R*8)  : W width.
C*       AME      : (R*8)  : electron mass.
C*       AMH      : (R*8)  : H mass.
C*       GMH      : (R*8)  : H width.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state nu spinor.
C*       SP2 (*)  : (C*16) : final-state nub spinor.
C*       SCH (*)  : (C*16) : final-state H wave function.
C* (Output)
C*       AMP(0)   : (C*16) : total t-channel amplitudes.
C*          (1)   : (C*16) : Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    97/02/09  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPNNH(GWF, GWWH, AMW, GMW, AME, AMH, GMH,
     .                         SPEM,SPEP, SP1,SP2, SCH, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
C--
      REAL   *8  GWF(2), GWWH, AMW, GMW, AME, AMH, GMH
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), SCH(3), AMP(0:1)
C--
      COMPLEX*16 WRK(6,2)
C--
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.D0)
         AMN   = 0
         GMN   = 0
         GME   = 0
       ENDIF
C--
C  Calculate internal lines.
C--
      CALL JIOXXX(SPEM,SP1,GWF,AMW,GMW, WRK(1, 1))
      CALL JIOXXX(SP2,SPEP,GWF,AMW,GMW, WRK(1, 2))
C--
C  Higgs.
C--
      CALL VVSXXX(WRK(1,2),WRK(1,1),SCH,GWWH, AMP(1))
C--
C  Sum up amplitudes.
C--
      AMP(0) = AMP(1)
C--
C  That's it.
C--
      RETURN
      END
