C        IMPLICIT REAL*8 ( A-H, O-Z )
C        COMMON /PAR/ YMAX, NY
C        AMT  = 150
C        RQCD = 0.365/0.2
C        R0   = 0.001/AMT
C        X    = LOG(RQCD/R0)
C        PRINT *, ' X = ', X
C  C--
C        YMAX = 10
C        NY   = 3000
C  C--
C        N = 100
C        XMIN = 0.001
C        XMAX = 10
C        DX   = (XMAX-XMIN)/N
C        WRITE(20,*) 'NEW FRAME'
C        WRITE(20,*) 'SET LIMIT X ', XMIN, ' ', XMAX
C        WRITE(20,*) 'SET SCALE X LOG'
C        WRITE(20,*) 'SET SCALE Y LOG'
C        WRITE(20,*) 'SET ORDER X Y'
C        DO 20 I = 0, N
C           X = XMIN + I*DX
C           F =  EI(-X)
C           WRITE(20,*) X, F
C  20    CONTINUE
C        WRITE(20,*) 'JOIN DOT'
C        DO 30 I = 0, N
C           X = XMIN + I*DX
C           F =  -EXP(-X)/X
C           F =  EI2(-X)
C           WRITE(20,*) X, F
C  30    CONTINUE
C        WRITE(20,*) 'JOIN'
C        END
C  C*
C  C*
C  C*
C
C        REAL FUNCTION EI2*8(X)
C
C        IMPLICIT REAL*8 ( A-H, O-Z )
C        REAL*8   X
C        DATA YMAX /30.D0/  NY /500/
C  C
C  C========< Entry Point >=============================================
C  C
C  C--
C  C  Calculate Ei(X).
C  C--
C        YMIN = -X
C        DY   = (YMAX-YMIN)/NY
C        SUM  = EXP( -YMIN )/YMIN + EXP( -YMAX )/YMAX
C        ALT  = 1
C        DO 10 IY = 1, NY - 1
C           Y   = YMIN + IY*DY
C           WT  = 3 + ALT
C           SUM = SUM + EXP( -Y )/Y*WT
C           ALT = -ALT
C  10    CONTINUE
C        EI2 = - SUM*DY/3
C  C--
C  C  That's it.
C  C--
C        RETURN
C        END
 
C*
C*  This routine calculates the exponential integral:
C*
C*                    _ /infty
C*                   /    dt
C*        Ei(x) = - /    --- exp(-t) ,
C*                _/-x    t
C*
C*  where x must be negative!
C*
 
      DOUBLE PRECISION FUNCTION EI(X)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   X
      DATA YMAX /10.D0/  NY /300/
C
C========< Entry Point >================================================
C
C--
C  Calculate Ei(X).
C--
      YMIN = LOG(-X)
      IF ( YMIN.GT.YMAX ) THEN
         EI = EXP(X)/X
      ELSE
         DY   = (YMAX-YMIN)/NY
         SUM  = EXP( -EXP(YMIN) ) + EXP( -EXP(YMAX) )
         ALT  = 1
         DO 10 IY = 1, NY - 1
            Y   = YMIN + IY*DY
            WT  = 3 + ALT
            SUM = SUM + EXP( -EXP(Y) )*WT
            ALT = -ALT
10       CONTINUE
         EI = - SUM*DY/3
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
