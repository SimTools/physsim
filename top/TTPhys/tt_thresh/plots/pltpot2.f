C*
C* This program draws 5-sigma bounds for alpha_s = 0.12 case.
C*
      IMPLICIT REAL*8 ( A-H, O-Z )
      COMMON /ABCD/   EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      REAL   *8       EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      COMMON /VPARA/  X(10)
      REAL   *8       X
C--
      REAL   *8       VDAT(0:6)
C--
      DATA NCALL /0/
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
C  Initialize parameters.
C--
      ALFS =  0.12D0
      AMT  =  150.D0
      AMB  =   5.D0
      AMZ  =  91.17D0
C--
C  Initialize input potential parameters.
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
C          = c_0 + C_1*log(r/r_o)*exp(-r/r_1) + a*r   ( r > r_0 )
C--
      SGR0   = 0.00954D0
      SGR1   = 0.35451D0
      SGSL   = 0.02467D0
C--
      R0     = 0.23530D0
      R1     = 0.37397D1
      SL     = 0.35719D0
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
C   Calculate lambda_4
C--
      CALL GTLAMB(BETA04,BETA14,A4,AMB,XLAM4)
C--
C  Prepare topdraw data.
C--
      NR   = 50
      RMN  = 0.01D0
      RMX  = 10.D0
      RLMN = LOG(RMN)
      RLMX = LOG(RMX)
      DRL  = (RLMX-RLMN)/NR
      DO 100 IR = 0, NR
         RL = RLMN + IR*DRL
         R  = EXP(RL)
         DO 10 I = 0, 6
C--
C  Modify potential parameters.
C--
            IF ( I.EQ.0 ) THEN
               DR0 = 0
               DR1 = 0
               DSL = 0
            ELSE IF ( I.EQ.1 ) THEN
               DR0 = -5*SGR0
               DR1 = 0
               DSL = 0
            ELSE IF ( I.EQ.2 ) THEN
               DR0 = +5*SGR0
               DR1 = 0
               DSL = 0
            ELSE IF ( I.EQ.3 ) THEN
               DR0 = 0
               DR1 = -5*SGR1
               DSL = 0
            ELSE IF ( I.EQ.4 ) THEN
               DR0 = 0
               DR1 = +5*SGR1
               DSL = 0
            ELSE IF ( I.EQ.5 ) THEN
               DR0 = 0
               DR1 = 0
               DSL = -5*SGSL
            ELSE IF ( I.EQ.6 ) THEN
               DR0 = 0
               DR1 = 0
               DSL = +5*SGSL
            ENDIF
            X(1)   = R0 + DR0
            X(2)   = R1 + DR1
            X(3)   = SL + DSL
            X(7)   = VP(X(1)) - X(3)*X(1)
            DLTR0  = X(1)**2*1.D-3
            VPP    = (VP(X(1))-VP(X(1)-DLTR0))/DLTR0
            X(8)   = X(1)*(VPP-X(3))/EXP(-X(1)/X(2))
C--
C  Calculate potential.
C--
            VDAT(I) = VRA(R)
10       CONTINUE
         WRITE(20,'(1E10.3,7F8.3)') R, (VDAT(I),I=0,6)
100   CONTINUE
C--
C  That's it.
C--
      END
