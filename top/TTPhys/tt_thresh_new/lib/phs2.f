CC**********************************************************************
C*  Real*8 function Phs2( q, p )
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION PHS2( Q, P )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       Q, P
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      REAL   *8       LAM, LAMAX, LAMIN, NTLAM2
C--
C  Statement function.
C--
      NTLAM2( LAM ) =
     .   0.0625D0*ETA**(-2)*LOG(2/(LAM*Y0))**2*Y0+0.0125D0*
     .   ETA**(-2)*LOG(2/(LAM*Y0))*Y0*(LAM**5*Y0*C1-2.5D0*
     .   LAM**4*C1+1.1111111D0*LAM**3*Y0*C2-3.3333333D0*LAM**2*
     .   C2+5*LAM*Y0-20)-(0.0022321429D0*ETA**(-2)*LAM*Y0
     .   )*(LAM**6*Y0**3*C1-9.3333333D0*LAM**5*Y0**2*C1+0.9333
     .   3333D0*LAM**4*Y0**3*C2+26.88D0*LAM**4*Y0*C1-9.3333333D0*
     .   LAM**3*Y0**2*C2-24.5D0*LAM**3*C1+2.3333333D0*LAM**2*Y0
     .   **3+29.037037D0*LAM**2*Y0*C2-28*LAM*Y0**2-28*LAM*
     .   C2+112*Y0)
C--
C
C===========< Entry Point >=============================================
C
C--
C  Calculate phase_2.
C--
      XSI = Q/GAMMAT
      ETA = P/GAMMAT
C--
      C1  = ( XSI**2 - ETA**2 )**2
      C2  = - 3*XSI**2 + ETA**2
      Y0  = 1/(1-RWT) * MT/GAMMAT * YCUT
C--
      LAMAX = 1/ABS( XSI - ETA )
      LAMIN = 1/( XSI + ETA )
C--
      IF  ( Y0*LAMIN .GT. 2.D0 )    THEN
            PHS2 = 0.
      ELSE IF ( Y0*LAMAX .GT. 2.D0 ) THEN
            PHS2 = NTLAM2( 2/Y0 ) - NTLAM2( LAMIN )
      ELSE
            PHS2 = NTLAM2( LAMAX ) - NTLAM2( LAMIN )
      END IF
C--
C  That's it.
C--
      RETURN
      END
