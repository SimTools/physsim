CC**********************************************************************
C*
C*=====================================================
C*  REAL*8 FUNCTION  ADJ2( EPS, STEP1, STEP2 )
C*=====================================================
C*
C* (Purpose)
C*    This function subroutine calculates the criterion adj2,
C*    which determines the upperbound of ! B_new/B_old - 1 ! .
C*    The criterion is adjusted so that the subroutine QCD1 would
C*    calculate Green's function for r<r_max, where r_max would be
C*    120% of r used in subroutine QCD2.
C* (Inputs)
C*       EPS   :REAL*8: criterion used for energy E.
C*       STEP1 :REAL*8: #steps required in subroutine QCD1 for E.
C*       STEP2 :REAL*8: #steps required in subroutine QCD2 for E.
C* (Output)
C*       ADJ2  :REAL*8: adjusted criterion to be used for E+delE.
C* (Update Record)
C*    92/04/14  K.Fujii                Original version by Sumino
C*                                     adopted to FACOM.
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION ADJ2( EPS, STEP1, STEP2 )
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL   *8  EPS
      INTEGER*4  STEP1, STEP2
C
C========< Entry Point >================================================
C
C--
C  Calculate adjusted criterion.
C--
      IF ( STEP1*10 .GT. STEP2*12 ) THEN
            ADJ2 = EPS/0.32D0
      ELSE
            ADJ2 = EPS*0.32D0
      END IF
C--
C  That's it.
C--
      RETURN
      END
