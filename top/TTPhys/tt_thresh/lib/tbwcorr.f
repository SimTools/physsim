CC*********************************************************************
C*========================
C*  REAL*8  FUNCTION H(Y)
C*========================
C*
C* (Purpose)
C*    This function subroutine calculates the factor necessary
C*    at the order alpha_s correction to the top decay width.
C*    The top width up to order alpha_s is given by
C*          Gammat = Gammat_zero*( 1 - CF*alpha_s/2/pi * h(r) )
C*    with  r = mW**2/mt**2.
C* (Input)
C*    Y     :REAL*8: (m_W/m_t)**2.
C* (Output)
C*    H     :REAL*8: O(alpha_s) correction to Gamma_t.
C* (Relation)
C*    Invokes DILOG.
C* (Update Record)
C*    93/04/14  K.Fujii                Original version by Sumino
C*                                     adopted to FACOM.
C*
CC*********************************************************************
 
      DOUBLE PRECISION FUNCTION H(Y)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL*8     Y
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.D0)
      ENDIF
C--
C  Calculate H(Y).
C--
      Z = 1 - Y
      H = PI**2 + 2*DILOG(Y) - 2*DILOG(Z)
     .          + ( 4*Y*(Z-2*Y**2)*LOG(Y)
     .               + 2*Z**2*(5+4*Y)*LOG(Z)
     .                  - Z*(5+9*Y-6*Y**2) )
     .                       /2/Z**2/(1+2*Y)
C--
C  That's it.
C--
      RETURN
      END
 
 
CC*********************************************************************
C*================================
C*  REAL*8  FUNCTION HBWG(R,YCUT)
C*================================
C*
C* (Purpose)
C*    This function subroutine calculates the factor necessary
C*    at the order alpha_s correction to the top quark
C*    partial width: t->bWg.
C*    The partial width at order alpha_s is given by
C*          Gamma(t->bWg) = Gammat_zero*CF*alpha_s/2/pi * h_bWg(r)
C*    with  r = mW**2/mt**2.
C* (Input)
C*    R     :REAL*8: (m_W/m_t)**2.
C*    YCUT  :REAL*8: cut on (m_bg/m_t)**2.
C* (Output)
C*    HBWG  :REAL*8: O(alpha_s) correction to Gamma(t->bWg).
C* (Relation)
C*    Invokes DILOG.
C* (Update Record)
C*    93/04/14  K.Fujii                Original version by Sumino
C*                                     adopted to FACOM.
C*
CC*********************************************************************
 
      DOUBLE PRECISION FUNCTION HBWG(R,YCUT)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL*8     R, YCUT
C
C========< Entry Point >================================================
C
C--
C  CALCULATE HBWG(R,YCUT)
C--
      Z = 1 - R
      B = - ( Z + YCUT )
      C = YCUT
      XP = ( - B + SQRT( B**2 - 4*C ) )/2
      XM = ( - B - SQRT( B**2 - 4*C ) )/2
C--
      HBWG = 2*( DILOG(XP) - DILOG(XM) - DILOG(XP/Z)
     .         + DILOG(XM/Z) - LOG(YCUT/Z**2)*LOG(XP/XM)/2 )
     .     + ( - 2*R*(1+R)*(1-2*R)*LOG( (1-XP)/(1-XM) )
     .         - (1+2*R)/2*( 7*Z**2 + 4*Z*YCUT + YCUT**2 )*LOG(XP/XM)
     .         + (XP-XM)*( 5 + 7*R - 8*R**2 + YCUT ) )/Z**2/(1+2*R)
C--
C  That's it.
C--
      RETURN
      END
 
 
C**********************************************************************
C   10/07/87 707161935  MEMBER NAME  LI2      (S)        M  FORTRAN77
C   10/07/87 707101837  MEMBER NAME  DILOG    (S)        M  FORTRAN77
C%]------------------- dilog.for
C%]
C-----------------------------------------------------------------------
C       Complex Spence Function  Li2(x)
C
C       If the argument is on the cut (x>1), It is assumed to lie above
C       the real axis (x+i0).
C
C       Ref.  G. 't Hooft and M. Veltman, Nucl. Phys. B153, 365 (1979).
C       Quoted accuracy is about 13 decimals in the worst case.
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION LI2(X)
      IMPLICIT COMPLEX*16(A-H,L,M,O-Z)
      REAL*8 DILOG,ZERO,ONE,HALF
      COMPLEX*16 IPI
      DATA IPI/(0D0,3.1415926535897932384D0)/
      DATA SP1/(1.644934066848226436D0,0D0)/
C                                  ] PI**2/6
      DATA ZERO/0D0/,ONE/1D0/,HALF/5D-1/
      IF(IMAG(X).EQ.ZERO.AND.ABS(REAL(X)).LE.ONE) THEN
        LI2=DILOG(REAL(X))
        RETURN
      ELSE IF(ABS(X).GT.ONE) THEN
        Y=1D0/X
        IND=1
        IF(IMAG(X).EQ.ZERO.AND.REAL(X).GT.ONE) IND=2
      ELSE
        Y=X
        IND=0
      END IF
      IF(REAL(Y).LE.HALF) THEN
        Z = - CLOG1(-Y)
        LI2 = CSP(Z)
      ELSE
        Z = - LOG(Y)
        LI2 = - CSP(Z) + SP1 + Z * CLOG1(-Y)
      END IF
      IF(IND.EQ.0) RETURN
      IF(IND.EQ.1) THEN
        LI2=-LI2 - HALF*LOG(-X)**2 - SP1
        RETURN
      ELSE
        LI2=-LI2 - HALF*LOG(X)**2 + IPI*LOG(X) + 2.*SP1
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C
      FUNCTION DILOG(X)
      IMPLICIT REAL*8 (A-H,L,M,O-Z)
      DATA SP1/1.644934066848226436D0/
C                                  ] PI**2/6
      IF(X.EQ.0.) THEN
        DILOG=.0
      ELSE IF(X.EQ.1.) THEN
        DILOG = SP1
      ELSE IF(X.EQ.-1.) THEN
        DILOG = -SP1*.5
      ELSE IF(ABS(X).LT.1.) THEN
        IF(X.LE..5) THEN
          Z = - LOG1(-X)
C         Z = - LOG(1.-X)
          DILOG = SP(Z)
        ELSE
          Z = - LOG(X)
          DILOG = - SP(Z) + SP1 + Z * LOG1(-X)
C         DILOG = - SP(Z) + SP1 + Z * LOG(1.-X)
        END IF
      ELSE
        DILOG=.0
        WRITE(6,*) 'ERROR in DILOG : X = ', X
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C
      FUNCTION CSP(Z)
      IMPLICIT REAL*8(A,B,D-H,M,O-Y)
      IMPLICIT COMPLEX*16(C,Z)
      DATA B0 /1D0/,
     1     B1 /-.25D0/,
     1     B2 / 2.7777777777777778D-02/,
     1     B4 /-2.7777777777777778D-04/,
     1     B6 / 4.7241118669690098D-06/,
     1     B8 /-9.1857730746619628D-08/,
     1     B10/ 1.8978869988970999D-09/,
     1     B12/-4.0647616451442253D-11/,
     1     B14/ 8.9216910204564524D-13/
C    1     B8 /-9.1857730746619636D-08/,
C    1     B10/ 1.8978869988970999D-09/,
C    1     B12/-4.0647616451442255D-11/,
C    1     B14/ 8.9216910204564526D-13/
C
      ZZ = Z * Z
      CSP = B0+Z*(B1+Z*(B2+ZZ*(B4+ZZ*(B6+ZZ*(B8+ZZ*(B10+ZZ*(B12+ZZ*B14)
     1     ))))))
      CSP = CSP * Z
C
      RETURN
      END
C----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION SP(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA B0 /1D0/,
     1     B1 /-.25D0/,
     1     B2 / 2.7777777777777778D-02/,
     1     B4 /-2.7777777777777778D-04/,
     1     B6 / 4.7241118669690098D-06/,
     1     B8 /-9.1857730746619628D-08/,
     1     B10/ 1.8978869988970999D-09/,
     1     B12/-4.0647616451442253D-11/,
     1     B14/ 8.9216910204564524D-13/
C    1     B8 /-9.1857730746619636D-08/,
C    1     B10/ 1.8978869988970999D-09/,
C    1     B12/-4.0647616451442255D-11/,
C    1     B14/ 8.9216910204564526D-13/
C
      XX = X * X
      SP = B0+X*(B1+X*(B2+XX*(B4+XX*(B6+XX*(B8+XX*(B10+XX*(B12+XX*B14)
     1     ))))))
      SP = SP * X
C
      RETURN
      END
C---------------------------------------------------------------------
      FUNCTION CLOG1(X)
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      COMPLEX*16 CLOG1,X
      DATA C1/1D0/,C2/.5D0/,C3/.33333333333333333D0/
      DATA C4/.25D0/,C5/.2D0/,C6/.166666666666666666D0/
      DATA C7/.142857142857142857D0/,C8/.125D0/
      DATA C9/.111111111111111111D0/,C10/.1D0/
      DATA C11/.09090909090909091D0/,C12/.08333333333333333333D0/
      IF(ABS(X).GT.1D-2) THEN
        CLOG1=LOG(1D0+X)
        RETURN
      ELSE
        CLOG1=X*(C1-X*(C2-X*(C3-X*(C4-X*(C5-X*(C6-X*(C7
     1       -X*(C8-X*(C9-X*(C10-X*(C11-X*C12)))))))))))
        RETURN
      ENDIF
      END
C-----------------------------------------------------------------------
C   14/07/87 707141333  MEMBER NAME  LOG1     (S)           FORTRAN77
      FUNCTION LOG1(X)
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DATA C1/1D0/,C2/.5D0/,C3/.33333333333333333D0/
      DATA C4/.25D0/,C5/.2D0/,C6/.166666666666666666D0/
      DATA C7/.142857142857142857D0/,C8/.125D0/
      DATA C9/.111111111111111111D0/,C10/.1D0/
      DATA C11/.09090909090909091D0/,C12/.08333333333333333333D0/
      IF(ABS(X).GT.1D-2) THEN
        LOG1=LOG(1D0+X)
        RETURN
      ELSE
        LOG1=X*(C1-X*(C2-X*(C3-X*(C4-X*(C5-X*(C6-X*(C7
     1       -X*(C8-X*(C9-X*(C10-X*(C11-X*C12)))))))))))
        RETURN
      ENDIF
      END
