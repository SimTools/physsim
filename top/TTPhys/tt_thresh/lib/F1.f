C        IMPLICIT REAL*8 ( A-H, O-Z )
C        N = 100
C        AMT = 150
C        AMW = 80
C        QCD = 0.2/0.365
C        EMIN = -20
C        EMAX = +20
C        DE   = (EMAX-EMIN)/N
C        WRITE(20,*) 'NEW FRAME'
C        WRITE(20,*) 'SET LIMIT X ', EMIN, ' ', EMAX
C        WRITE(20,*) 'SET ORDER X Y'
C        DO 10 I = 0, N
C           E = EMIN + I*DE
C           QCUT = SQRT( (3*AMT+E-AMW)*(3*AMT+E+AMW)
C       .            *(AMT+E-AMW)*(AMT+E+AMW) )/2/(2*AMT+E)
C           WRITE(20,*) E, QCUT/QCD
C  10    CONTINUE
C        WRITE(20,*) 'JOIN'
C  C--
C        XMIN = 10
C        XMAX = 1510
C        DX   = (XMAX-XMIN)/N
C        WRITE(20,*) 'NEW FRAME'
C        WRITE(20,*) 'SET LIMIT X ', XMIN, ' ', XMAX
C        WRITE(20,*) 'SET ORDER X Y'
C        DO 20 I = 0, N
C           X = XMIN + I*DX
C           F =  F1(X)
C           WRITE(20,*) X, F
C  20    CONTINUE
C        WRITE(20,*) 'JOIN'
C        END
C*
C*
C*
 
      DOUBLE PRECISION FUNCTION F1(X)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   X
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI     = ACOS(-1.D0)
      END IF
C--
C  Calculate F_x(x).
C--
      NY   = 3000
      YMAX = 3
      DY   = YMAX/NY
      SUM  = 0
      DO 100 IY = 0, NY
         Y   = IY*DY
         IF ( IY.EQ.0 .OR. IY.EQ.NY ) THEN
            WT = 1
         ELSE IF ( MOD(IY,2).EQ.0 ) THEN
            WT = 2
         ELSE
            WT = 4
         ENDIF
         EXPY = EXP(Y)
         SUM = SUM + Y*EXPY*SIN( X*EXP( -EXPY ) )*WT
100   CONTINUE
      F1  = 2/PI*SUM*DY/3
C>>>
C     PRINT *, ' X F1 = ', X, F1
C>>>
C--
C  That's it.
C--
      RETURN
      END
