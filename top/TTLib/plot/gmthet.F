      SUBROUTINE GMTHET(Q,GM)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       Q, GM
C--
      COMMON /KFFLAG/ IPOT, IRGM
      INTEGER*4       IPOT, IRGM
      COMMON /PARAM1/ E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      REAL   *8       E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
C--
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
      IF (NCALL.EQ.0) THEN
         NCALL = 1
         PI  = ACOS(-1.D0)
         SQH = SQRT(0.5D0)
      ENDIF
C--
C  Calculate Gamma_theta
C--
      GF  = (PI/SQRT(2.D0))*(ALF/S2W)/AMW**2
      W   = AMW/AMT
      C10 = 1 - 3*W**4 + 2*W**6
      IF (IRGM.EQ.1) THEN
         C1E = 0
         C1  = C10 + C1E
         C2  = C10/2
#if 1
C        non-relativistic time retardation
C        gamma^{-1} ~ 1 - (Q/m)**2/2
         RET = C1 - C2*Q**2/AMT/AMT
#else
C        relativistic time retardation
         e1  = sqrt(q**2 + amt**2)
         ret = C10*amt/e1
#endif
      ELSE
         RET = C10
      ENDIF
      GM0 = VTB**2*GF*SQH*(AMT**3/4/PI)*C10
      GM  = VTB**2*GF*SQH*(AMT**3/4/PI)*RET
#if 0
      print *, 'GMTHET = ', gm, ' GMTHET0 = ', gm0
#endif
C--
C  That's it.
C--
      RETURN
      END
