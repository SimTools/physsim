CC**********************************************************************
C*
C*=====================================================
C*  REAL*8 FUNCTION ADJ1( DROLD, STPOLD, DR, STEP )
C*=====================================================
C*
C* (Purpose)
C*    This function subroutine calculates the adjusted integration
C*    step in solving Schroedinger equation, so that the produced
C*    Green's function has about 250 (=ideal) points.
C* (Inputs)
C*       DROLD  :REAL*8: integration step used for energy E-delE.
C*       STPOLD :REAL*8: number of steps required for solving
C*                       Green's function for energy E-delE.
C*       DR     :REAL*8: integration step used for energy E.
C*       STEP   :REAL*8: number of steps required for solving
C*                         Green's function for energy E.
C* (Output)
C*       ADJ1   :REAL*8: adjusted integration step to be used for
C*                       energy E+delE.
C* (Update Record)
C*    92/04/14  K.Fujii                Original version by Sumino
C*                                     adopted to FACOM.
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION ADJ1( DROLD, STPOLD, DR, STEP )
 
      REAL   *8  DROLD, DR, RMAX
      INTEGER*4  STPOLD, STEP, IDEAL
      PARAMETER  ( IDEAL = 250 )
C
C========< Entry Point >================================================
C
C--
C  Calculate adjusted step size.
C--
      RMAX = 2*DR*DBLE(STEP) - DROLD*DBLE(STPOLD)
      ADJ1 = RMAX/DBLE(IDEAL)
C--
C  That's it.
C--
      RETURN
      END
