C
C        IMPLICIT  REAL*8 ( A-H, O-Z )
C        REAL*8    A(3)
C        DATA  A   / 5.D0, -6.D0, 1.D0 /
C        DATA NPAR / 3 /
C
C        DO 10 I = 1, 10
C           PRINT *, ' Input X0, Y0 '
C           READ(5,*) X0, Y0
C           CALL USOLVE(NPAR,A,X0,Y0,X,Y,IRT)
C           PRINT *, ' X, Y = ', X, Y
C  10    CONTINUE
C        END
 
      SUBROUTINE USOLVE(NPAR,A,X04,Y04,X4,Y4,IRT)
 
      IMPLICIT  REAL*8 ( A-H, O-Z )
      INTEGER*4    NPAR, IRT
      REAL*4       A(0:NPAR-1), X04, Y04, X4, Y4
      DATA   TEST / 1.D-6 /  NTRYMX  / 10000 /
C
C========< Entry Point >================================================
C
C--
C  Reset return flag.
C--
      IRT    = 0
C--
C  Set initial values.
C--
      EPS    = 1.D-3
      ADYS   = 1.D20
      X0     = X04
      Y0     = Y04
      X      = X0
C--
C  Start Newton's method.
C--
      NTRY = 0
1     NTRY = NTRY + 1
C--
      Y    = 0
      DO 20 IE = NPAR-1, 0, -1
         Y = Y*X + A(IE)
20    CONTINUE
C--
      DYDX = 0
      DO 30 IE = NPAR-1, 1, -1
         DYDX = DYDX*X + IE*A(IE)
30    CONTINUE
C--
      ADYS  = ADY
      DY    = Y - Y0
      ADY   = ABS(DY)
      IF ( ADY.GT.ADYS ) THEN
         EPS = EPS*1.D+1
      ELSE
         EPS = EPS*1.D-1
      ENDIF
      DX     = -DY/DYDX/(1.+EPS)
      IF ( ADY.GT.TEST ) THEN
         X   = X + DX
         IF ( NTRY.LT.NTRYMX )                   GO TO 1
         PRINT *, ' USOLVE iteration over'
         IRT = NTRY
      ENDIF
      X4 = X
      Y4 = Y
      RETURN
      END
