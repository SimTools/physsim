C*
C* (Update Record)
C*   92/03/19  K.Fujii       Parametrization of Cho's new potential
C*                           parameters.
C*   92/03/19  K.Fujii       This version uses the exact values
C*                           obtained by Cho.
C*
 
      SUBROUTINE GTPPAR(X,PPAR)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   X, PPAR(3)
C--
      PARAMETER      ( NPT = 5 )
      REAL*4           YDAT(NPT,3), DYDAT(NPT,3)
      DATA YDAT  / 0.20823, 0.23347, 0.23530, 0.20681, 0.16754,
     .             3.9620 , 3.8077 , 3.7397 , 3.7453 , 3.7485 ,
     .             0.35914, 0.35398, 0.35719, 0.37349, 0.40100 /
      DATA DYDAT / 0.01028, 0.01095, 0.00954, 0.00572, 0.00203,
     .             0.27724, 0.28706, 0.35451, 0.40334, 0.29333,
     .             0.02136, 0.02130, 0.02466, 0.02769, 0.02078 /
C
C==================< Entry Point >======================================
C
C--
C  Calculate PPAR.
C--
      IF ( X.GE.0.099D0 .AND. X.LE.0.101D0 ) THEN
         J = 1
      ELSE IF ( X.GE.0.109D0 .AND. X.LE.0.111D0 ) THEN
         J = 2
      ELSE IF ( X.GE.0.119D0 .AND. X.LE.0.121D0 ) THEN
         J = 3
      ELSE IF ( X.GE.0.129D0 .AND. X.LE.0.131D0 ) THEN
         J = 4
      ELSE IF ( X.GE.0.139D0 .AND. X.LE.0.141D0 ) THEN
         J = 5
      ELSE
         PRINT *, ' ERROR in GTPPAR: ALFS = ', X
      ENDIF
C>>>
C     PPAR(1) = YDAT(J,1) - 5*DYDAT(J,1)
C     PPAR(1) = YDAT(J,1)
      PPAR(1) = YDAT(J,1) + 5*DYDAT(J,1)
C>>>
      PPAR(2) = YDAT(J,2)
      PPAR(3) = YDAT(J,3)
C--
C  That's it.
C--
      RETURN
      END
