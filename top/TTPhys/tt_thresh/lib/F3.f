C
C        IMPLICIT REAL*8 ( A-H, O-Z )
C  C--
C        N    = 100
C        XMIN = 0.1
C        XMAX = 100.1
C        DX   = (XMAX-XMIN)/N
C        WRITE(20,*) 'NEW FRAME'
C        WRITE(20,*) 'SET LIMIT X ', XMIN, ' ', XMAX
C        WRITE(20,*) 'SET ORDER X Y'
C        DO 20 I = 0, N
C           X = XMIN + I*DX
C           F =  F3(X)
C           WRITE(20,*) X, F
C  20    CONTINUE
C        WRITE(20,*) 'JOIN'
C        END
C*
C*
C*
 
      DOUBLE PRECISION FUNCTION F3(X)
 
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
      YMAX = 50
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
         SUM = SUM + Y*SIN( X/EXPY )*WT
100   CONTINUE
      F3  = 2/PI*SUM*DY/3
C>>>
C     PRINT *, ' X F3 = ', X, F3
C>>>
C--
C  That's it.
C--
      RETURN
      END
