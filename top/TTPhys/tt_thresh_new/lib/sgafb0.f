 
      SUBROUTINE SGAFB0(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,
     .                  E,SIGTOT,DLTAFB)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
      INTEGER* 4        MODE
      REAL*8            ALPS, ALP, SN2W, AMSZ, AMSW, AMSB, AMST, VTB2,
     .                  E, SIGTOT, DLTAFB
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      CALL GRQINT(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,E)
C--
C  Calculate Green's functions.
C--
      CALL GRQSGA(MODE,E,SIGTOT,DLTAFB)
C--
C  That's it.
C--
      RETURN
      END
