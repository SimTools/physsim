C  *INCLUDE F3
C        IMPLICIT REAL*8 ( A-H, O-Z )
C  C--
C        NX   = 100
C        XMIN =   0
C        XMAX =  60
C        DX   = (XMAX-XMIN)/NX
C  C--
C        WRITE(20,*) 'NEW FRAME'
C        WRITE(20,*) 'SET LIMIT X ', XMIN, ' ', XMAX
C        WRITE(20,*) 'SET LIMIT Y 0 7'
C        WRITE(20,*) 'SET ORDER X Y'
C        DO 10 IX = 0, NX
C           X  = XMIN + IX*DX
C           Y  = F3(X)
C           WRITE(20,*) X, Y
C  10    CONTINUE
C        WRITE(20,*) 'JOIN DOT'
C  C--
C        NX   =  300
C        DX   = (XMAX-XMIN)/NX
C        DO 20 IX = 0, NX
C           X  = XMIN + IX*DX
C           Y  = FH(X)
C           WRITE(20,*) X, Y
C  20    CONTINUE
C        WRITE(20,*) 'JOIN'
C        END
C*
C*
C*
 
      DOUBLE PRECISION FUNCTION FH(X)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   X
C--
      REAL*8 FHDATA(0:60)
      DATA NX     /        60 /
      DATA XMIN   /.00000D+00 /
      DATA DX     /.10000D+01 /
      DATA FHDATA /
     . 0.0000D+00, 0.6250D+00, 0.1185D+01, 0.1638D+01, 0.1972D+01,
     . 0.2208D+01, 0.2379D+01, 0.2519D+01, 0.2648D+01, 0.2770D+01,
     . 0.2882D+01, 0.2980D+01, 0.3065D+01, 0.3141D+01, 0.3213D+01,
     . 0.3283D+01, 0.3350D+01, 0.3412D+01, 0.3469D+01, 0.3522D+01,
     . 0.3572D+01, 0.3620D+01, 0.3668D+01, 0.3714D+01, 0.3756D+01,
     . 0.3796D+01, 0.3835D+01, 0.3872D+01, 0.3909D+01, 0.3945D+01,
     . 0.3979D+01, 0.4012D+01, 0.4043D+01, 0.4073D+01, 0.4103D+01,
     . 0.4133D+01, 0.4161D+01, 0.4188D+01, 0.4215D+01, 0.4240D+01,
     . 0.4266D+01, 0.4291D+01, 0.4315D+01, 0.4339D+01, 0.4361D+01,
     . 0.4384D+01, 0.4406D+01, 0.4427D+01, 0.4449D+01, 0.4469D+01,
     . 0.4489D+01, 0.4509D+01, 0.4528D+01, 0.4547D+01, 0.4566D+01,
     . 0.4585D+01, 0.4603D+01, 0.4620D+01, 0.4637D+01, 0.4655D+01,
     . 0.4672D+01/
C
C========< Entry Point >================================================
C
C--
C  Calculate f_x(x) = f_1(x) + b_1/b_0^2*f_2(x).
C--
      IX = (X-XMIN)/DX
      IF ( IX.LT.0 .OR. IX.GE.NX ) THEN
         PRINT *, ' Error in FH: X = ', X, ' out of range.'
      ELSE
         FH = FHDATA(IX) + (FHDATA(IX+1)-FHDATA(IX))*(X-XMIN-IX*DX)/DX
      END IF
C--
C  That's it.
C--
      RETURN
      END
