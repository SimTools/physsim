CC**********************************************************************
C*  Real*8 function PHS0( q, p )
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION PHS0( Q, P )
 
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
      REAL   *8       LAM, LAMAX, LAMIN, NTLAM0
C--
C  Statement funtion.
C--
      NTLAM0( LAM ) =
     .     0.027777778D0*LOG(2/(LAM*Y0))*LAM**2*Y0*(KAPPA**2
     .     *LAM*Y0-3*KAPPA**2+9)-(0.0041666667D0*LAM**2*Y0)*
     .     (KAPPA**2*LAM**3*Y0**3-10*KAPPA**2*LAM**2*Y0**2+
     .     31.111111D0*KAPPA**2*LAM*Y0-30*KAPPA**2-20*LAM*Y0+
     .     30)
C--
C
C===========< Entry Point >=============================================
C
C--
C  Calculate phase_0.
C--
      XSI   = Q/GAMMAT
      ETA   = P/GAMMAT
C--
      Y0    = 1/(1-RWT) * MT/GAMMAT * YCUT
      LAMAX = 1/ABS( XSI - ETA )
      LAMIN = 1/( XSI + ETA )
C--
      IF  ( Y0*LAMIN .GT. 2.D0 )    THEN
            PHS0 = 0
      ELSE IF ( Y0*LAMAX .GT. 2.D0 ) THEN
            PHS0 = NTLAM0(  2/Y0 ) - NTLAM0( LAMIN )
      ELSE
            PHS0 = NTLAM0( LAMAX ) - NTLAM0( LAMIN )
      END IF
C--
C  That's it.
C--
      RETURN
      END
