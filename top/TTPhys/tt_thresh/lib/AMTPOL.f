      REAL*8 FUNCTION AMTPOL(MODE, AMTMSB, ALFSZ, AMZ )    
   
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4 MODE
      REAL*8    AMTMSB, ALFSZ, AMZ
      DATA      NCALL / 0 /
C  
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL .EQ. 0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.D0)
      ENDIF
C--
C  Calculate pole mass.
C--
      IF (MODE.EQ.1) THEN
         ALRM   = LOG(AMZ/AMTMSB)
         AMTPOL = AMTMSB * ( 1 + ALFSZ * ( 4/(3*PI)
     .                      + ALFSZ * ( 0.834515 + 0.517864*ALRM
     .                      + ALFSZ * ( 2.36827  + 2.24437*ALRM
     .                                           + 0.631891*ALRM**2 ))))
      ELSE
         AMTPOL = AMTMSB
      ENDIF
C--
C  That's it.
C--
      RETURN
      END         


