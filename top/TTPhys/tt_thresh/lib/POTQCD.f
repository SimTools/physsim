CC**********************************************************************
C*
C*========================================----------------===
C* Subroutine SETVPR(ALFS,AMT,AMB,AMZ,GMT,ALFST,ALFSP,RQCD)
C*========================================----------------===
C*
C* (Function)
C*    Initialize potential parameters and tabulate QCD potential for
C*    later use in V(R) calculation.
C*     -Usage-
C*          CALL SETVPR(...)
C*          VR = V(R)
C* (Inputs)
C*    ALFS  :(R*8): alpha_s(m_Z).
C*    AMT   :(R*8): m_t.
C*    AMB   :(R*8): m_b.
C*    AMZ   :(R*8): m_Z.
C*    GMT   :(R*8): gamma_t.
C* (Outputs)
C*    ALFST :(R*8): alpha_s(m_t).
C*    ALFSP :(R*8): alpha_s(m_t*alpha_s(m_Z)).
C*    RQCD  :(R*8): r_QCD.
C* (Update Record)
C*    91/10/23  K.Fujii       Derived from PARAV1 by Cho.
C*                            Modified to use the parametrization of
C*                            potential parameters to allow arbitrary
C*                            alpha_s values. GMT is used only for
C*                            determining the range of tabulation.
C*    92/03/18  K.Fujii       New potential parametrization.
C*
CC**********************************************************************
 
      SUBROUTINE SETVPR(ALFS,AMT,AMB,AMZ,GMT,ALFST,ALFSP,RQCD)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       ALFS, AMT, AMB, AMZ, GMT, ALFST, ALFSP, RQCD
C--
      COMMON /ABCD/   EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      REAL   *8       EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      COMMON /VPARA/  X(10)
      REAL   *8       X
C--
      PARAMETER      ( NTRMAX = 3000 )
C     PARAMETER      ( NTRMAX = 1000 )
      COMMON /VR0/    VR(0:NTRMAX), TR(0:NTRMAX), DTR
      REAL   *8       VR          , TR          , DTR
C--
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Calculate QCD parameters.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI     = ACOS(-1.D0)
         EUL    = 0.5772D0
         BETA04 = (33.D0-2.D0*4.D0)/6
         BETA14 = (153.D0-19.D0*4.D0)/12
         BETA05 = (33.D0-2.D0*5.D0)/6
         BETA15 = (153.D0-19.D0*5.D0)/12
         BLER4  = (BETA04+1.D0/3)*EUL+93.D0/37-5.D0/6
         BLER5  = (BETA05+2.D0/3)*EUL+93.D0/37-5.D0/6
         A      = BETA05*EUL + 93.D0/37 - 25.D0/18
      END IF
C--
C  Calculate input potential parameters.
C     X(1) : r_0
C      (2) : r_1
C      (3) : a
C      (4) : m_c
C      (5) : m_b
C      (6) : m_t
C      (7) : c_0
C      (8) : c_1
C  where these determines potential as follows:
C     V(r) = V_p(r)                                   ( r < r_0 )
C          = c_0 + c_1*log(r/r_o)*exp(-r/r_1) + a*r   ( r > r_0 )
C--
      CALL GTPPAR(ALFS,X(1))
      X(4)   = 1.5D0
      X(5)   = AMB
      X(6)   = AMT
C--
C  Calculate lambda_5 from alphaz.
C--
      A5    = ALFS/PI
      CALL GTLAMB(BETA05,BETA15,A5,AMZ,XLAM5)
C--
C  Calculate alphas (nf=5) at mu=mb
C--
      CALL GTALFS(BETA05,BETA15,XLAM5,AMB,A4)
C--
C  Calculate lambda_4.
C--
      CALL GTLAMB(BETA04,BETA14,A4,AMB,XLAM4)
C--
C  Calculate alphas (nf=5) at mu=mt
C--
      CALL GTALFS(BETA05,BETA15,XLAM5,AMT,A5)
      ALFST = A5*PI
C--
C  Calculate alphas (nf=5) at mu=mt*alpha_s(m_Z).
C--
      CALL GTALFS(BETA05,BETA15,XLAM5,AMT*ALFS,A5P)
      ALFSP = A5P*PI
C--
C  Calculate r_QCD.
C--
      RQCD = 1/XLAM5 * EXP(-A/BETA05)
     .               * ( 2*BETA15/BETA05**2 )**(BETA15/BETA05**2)
C--
C  Now calculate c_0 and c_1.
C--
      R0     = X(1)
      X(7)   = VP(R0) - X(3)*R0
      DR0    = R0**2*1.D-3
      VPP    = (VP(R0)-VP(R0-DR0))/DR0
      X(8)   = R0*(VPP-X(3))/EXP(-R0/X(2))
C--
C  Tablulate QCD potential.
C--
C>>>
      PRINT *, 'SETVPR now tabulates QCD potential.'
C>>>
C     RMIN = .00005/AMT
      RMIN = .00001/AMT
C>>>
      IF ( GMT.EQ.0.D0 ) THEN
         RMAX = 1.D0
      ELSE
         RMAX = 30.D0/SQRT(AMT*GMT)*2.
      END IF
      TRMIN = LOG(RMIN)
      TRMAX = LOG(RMAX)
      DTR   = (TRMAX-TRMIN)/NTRMAX
      DO 300 ITR = 0, NTRMAX
         TR(ITR) = TRMIN + ITR*DTR
         RR      = EXP(TR(ITR))
         VR(ITR) = VRA(RR)
CX       WRITE(10,*) RR, VR(ITR)
300   CONTINUE
C>>>
      PRINT *, 'Tabulation ended.'
C>>>
C--
C  That's it.
C--
      RETURN
      END
 
CC**********************************************************************
C*
C*===============================--===
C* Subroutine GTLAMB(BT0,BT1,A,Q,AL)
C*===============================--===
C*
C* (Function)          __
C*    Calculate lambda_MS(2) from a = alpha_s(Q)/pi.
C* (Inputs)
C*    BT0   :(R*8): beta_0.
C*    BT1   :(R*8): beta_1.
C*    A     :(R*8): alpha_s(Q)/pi.
C*    Q     :(R*8): energy scale.
C* (Output)                __
C*    AL    :(R*8): lambda_MS(2)
C* (Update Record)
C*    91/10/29  K.Fujii       Original version.
C*
CC**********************************************************************
 
      SUBROUTINE GTLAMB(BT0,BT1,A,Q,AL)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   BT0, BT1, A, Q, AL
C
C========< Entry Point >================================================
C
C--
C  Calculate lambda from a = alpha_s(Q)/pi.
C--
      FF1  = -1/BT0/A
      FF2  = BT1/BT0**2*LOG(2/BT0*(1/A+BT1/BT0))
      AL   = Q*EXP(FF1+FF2)
C--
C  That's it.
C--
      RETURN
      END
 
CC**********************************************************************
C*
C*===============================--===
C* Subroutine GTALFS(BT0,BT1,A,Q,AL)
C*===============================--===
C*
C* (Function)                                      __
C*    Calculate a = alpha_s(Q)/pi for given lambda_MS(2).
C* (Inputs)
C*    BT0   :(R*8): beta_0.
C*    BT1   :(R*8): beta_1.__
C*    AL    :(R*8): lambda_MS(2)
C*    Q     :(R*8): energy scale.
C* (Output)
C*    A     :(R*8): alpha_s(Q)/pi.
C* (Update Record)
C*    91/10/29  K.Fujii       Original version.
C*
CC**********************************************************************
 
      SUBROUTINE GTALFS(BT0,BT1,AL,Q,A)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   BT0, BT1, AL, Q, A
      REAL*8   WF(150)
C
C========< Entry Point >================================================
C
C--
C  Calculate alphas at Q.
C--
      GG2   = BT0**2/BT1
      GG4   = (GG2/2)**(1/GG2)
      GG5   = BT0**2*LOG(Q/AL)
      GG8   = BT1*LOG(2/GG2)
      PP    = BT1/(GG5+GG8)
      WF(1) = PP
      P4    = 0
      DO 100 JK = 1, 149
         WF(JK+1) = PP/(1+PP*LOG((WF(JK)+1)/WF(JK)))
         IF ( ABS(WF(JK+1)-WF(JK)).LT.1.0D-04 ) THEN
            P4 = WF(JK+1)
                                                 GO TO 150
         END IF
100   CONTINUE
150   IF ( P4.EQ.0.D0 ) WRITE(6,*) 'ALPHAS_4 ERROR'
      A    = BT0*P4/BT1
C--
C  That's it.
C--
      RETURN
      END
 
C*
C*   This routine interpolates the potential values tabulated in
C*   SETVPR by calling VRA.
C*
 
      DOUBLE PRECISION FUNCTION V(R)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       R
C--
      PARAMETER      ( NTRMAX = 3000 )
C     PARAMETER      ( NTRMAX = 1000 )
      COMMON /VR0/    VR(0:NTRMAX), TR(0:NTRMAX), DTR
      REAL   *8       VR          , TR          , DTR
C
C========< Entry Point >================================================
C
C--
C  Interpolate potential table.
C--
      TRMIN = TR(0)
      T     = DLOG(R)
      ITR   = (T-TRMIN)/DTR
      A1    = T - TR(ITR)
      A2    = TR(ITR+1) - T
      IF ( DABS(A1).LT.1.D-12 ) A1 = 0.D0
      IF ( DABS(A2).LT.1.D-12 ) A2 = 0.D0
      IF ( A1*A2.LT.0.D0 ) THEN
         PRINT *, ' Error in V(r): ', ITR, R, A1, A2
         PRINT *, '          R   =', R
         PRINT *, '          RMIN =', EXP(TRMIN)
      END IF
      V     = (A2*VR(ITR)+A1*VR(ITR+1))/DTR
C--
C  That's it.
C--
      RETURN
      END
 
C*
C*   This function calculates the 2-loop QCD potential as
C*           VRA(r) = V_p(r)  : r < r_0
C*                  = V_IL(r) : r > r_0
C*
 
      DOUBLE PRECISION FUNCTION VRA(R)
 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL   *8       R
      COMMON /VPARA/  X(10)
      REAL   *8       X
C
C========< Entry Point >================================================
C
C--
C  Store potential parameters to local varialbles.
C--
      R0 = X(1)
      R1 = X(2)
      A  = X(3)
      C0 = X(7)
      C1 = X(8)
C--
C  Branch on the r value.
C--
      IF ( R.LE.R0 ) THEN
C--
C  Calculate the short distance potential.
C--
         VRA = VP(R)
C--
C  Calculate intermediate and long distance potential.
C--
      ELSE
         VRA  = C0 + C1*LOG(R/R0)*EXP(-R/R1) + A*R
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
 
C*
C*   This function calculates the 2-loop QCD short distance
C*   potential.
C*
 
      DOUBLE PRECISION FUNCTION VP(R)
 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL   *8       R
C--
      COMMON /ABCD/   EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      REAL   *8       EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      COMMON /VPARA/  X(10)
      REAL   *8       X
C--
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Calculate numerical constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.D0)
         EX56  = EXP( 5.D0/6 )
      ENDIF
C--
C  Set charm and bottom masses.
C--
      XMC   = X(4)
      XMB   = X(5)
C--
      XMCR  = XMC*R
      XMBR  = XMB*R
C--
C  Calculate the exponential integrals.
C--
      EIC   = EI( -EX56*XMCR )
      EIB   = EI( -EX56*XMBR )
C--
C  Calculate A(r) and determine the scale Q.
C--
      IF ( XMBR.LT.1.D-2 ) THEN
         AR = BLER5 - 2*EUL/3 - 5.D0/9
      ELSE IF ( XMCR.LT.1.D-2 ) THEN
         AR = BLER5 - EUL/3 - 5.D0/18 + ( LOG(XMBR) - EIB )/3
      ELSE
         AR = BLER5 + ( LOG(XMCR) - EIC + LOG(XMBR) - EIB  )/3
      END IF
      Q     = EXP(-AR/BETA05)/R
C--
C  Calculate the short distance potential.
C--
      CALL GTALFS(BETA05,BETA15,XLAM5,Q,A5)
      ALFS  = PI*A5
      VP    = -4*ALFS/(3*R)
C--
C  That's it.
C--
      RETURN
      END
 
C*
C* (Update Record)
C*   92/03/19  K.Fujii       Parametrization of Cho's new potential
C*                           parameters.
C*
 
      SUBROUTINE GTPPAR(X,PPAR)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   X, PPAR(3)
      REAL*8   A(3,3)
      DATA (A(I,1),I=1,3)
     .    /-0.11540E+01, 0.24173E+02,-0.10510E+03/
      DATA (A(I,2),I=1,3)
     .    / 0.84780E+01,-0.73986E+02, 0.28710E+03/
      DATA (A(I,3),I=1,3)
     .    / 0.10357E+01,-0.12334E+02, 0.55747E+02/
C
C==================< Entry Point >======================================
C
C--
C  Calculate PPAR.
C     PPAR(1) = r_0
C         (2) = r_1
C         (3) = a
C--
      DO 100 J = 1, 3
         Y = A(3,J)
         DO 10 I = 2, 1, -1
            Y = Y*X + A(I,J)
10       CONTINUE
         PPAR(J) = Y
100   CONTINUE
C--
C  That's it.
C--
      RETURN
      END
