CINCLUDE EI
C*INCLUDE QCDPOT1
CINCLUDE POTQCD
C*
C* This program draws alpha_s dependence.
C*
      IMPLICIT REAL*8 ( A-H, O-Z )
C--
      PARAMETER      ( NPT = 5 )
      REAL*4          XDAT(NPT), YDAT(NPT,3), DYDAT(NPT,3)
C--
      COMMON /ABCD/   EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
      REAL   *8       EUL, XLAM4, XLAM5, BETA04, BETA14, BETA05, BETA15,
     .                BLER4, BLER5
C--
      COMMON /VPARA/  X(10)
      REAL   *8       X
C--
      DATA XDAT  / 0.10   , 0.11,    0.12,    0.13,    0.14    /
      DATA YDAT  / 0.20823, 0.23347, 0.23530, 0.20681, 0.16754,
     .             3.9620 , 3.8077 , 3.7397 , 3.7453 , 3.7485 ,
     .             0.35914, 0.35398, 0.35719, 0.37349, 0.40100 /
      DATA DYDAT / 0.01028, 0.01095, 0.00954, 0.00572, 0.00203,
     .             0.27724, 0.28706, 0.35451, 0.40334, 0.29333,
     .             0.02136, 0.02130, 0.02466, 0.02769, 0.02078 /
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
C  Loop over alpha_s.
C--
C>>>
C     DO 1000 IALFS = 1, 4
      DO 1000 IALFS = 3, 3
C>>>
         ALFS = XDAT(IALFS)
         WRITE(20,'(''( alpha_s = '',F5.2)') ALFS
         WRITE(20,'(''SET ORDER X Y'')')
C--
C  Initialize parameters.
C--
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
         X(1)   = YDAT(IALFS,1)
         X(2)   = YDAT(IALFS,2)
         X(3)   = YDAT(IALFS,3)
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
C  Calculate lambda_4
C--
         CALL GTLAMB(BETA04,BETA14,A4,AMB,XLAM4)
C--
C  Calculate remaining parameters.
C--
         X(7)   = VP(X(1)) - X(3)*X(1)
         DLTR0  = X(1)**2*1.D-3
         VPP    = (VP(X(1))-VP(X(1)-DLTR0))/DLTR0
         X(8)   = X(1)*(VPP-X(3))/EXP(-X(1)/X(2))
C--
C  Prepare topdraw data.
C--
C>>>
C        NR   = 50
C        RMN  = 0.01D0
         NR   = 100
         RMN  = 0.001D0
C>>>
         RMX  = 10.D0
         RLMN = LOG(RMN)
         RLMX = LOG(RMX)
         DRL  = (RLMX-RLMN)/NR
         DO 100 IR = 0, NR
            RL = RLMN + IR*DRL
            R  = EXP(RL)
C--
C  Calculate potential.
C--
            VDT    = VRA(R)
            WRITE(20,'(2E15.5)') R, VDT
100      CONTINUE
         WRITE(20,'(''JOIN'')')
1000  CONTINUE
C--
C  That's it.
C--
      END
