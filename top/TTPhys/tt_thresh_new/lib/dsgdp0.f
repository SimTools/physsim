 
      SUBROUTINE DSGDP0(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,
     .                  E, Q,DSG)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
      INTEGER* 4        MODE
      REAL*8            ALPS, ALP, SN2W, AMSZ, AMSW, AMSB, AMST, VTB2,
     .                  E, Q, DSG
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      CALL GRQINT(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,E)
C--
C  Calculate differential cross section.
C--
      CALL GRQDSG(MODE,E,Q,DSG)
C--
C  That's it.
C--
      RETURN
      END
