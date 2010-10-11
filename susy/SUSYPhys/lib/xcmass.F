 
      SUBROUTINE XCMASS(AM2,AMU,BT,AMW,AMXC,PHIL,PHIR,EPSR)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     AM2, AMU, BT, AMW, AMXC(2), PHIL, PHIR, EPSR
      REAL*8     AMSMAT(2,2), UL(2,2), UR(2,2)
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.)
         SQ2   = SQRT(2.)
      ENDIF
C--
C  Calculate mass matrix.
C--
      CSB  = COS(BT)
      SNB  = SIN(BT)
      SN2B = 2*SNB*CSB
      CS2B = (CSB-SNB)*(CSB+SNB)
C--
      AMSMAT(1,1) = AM2
      AMSMAT(1,2) = SQ2*AMW*CSB
      AMSMAT(2,1) = SQ2*AMW*SNB
      AMSMAT(2,2) = AMU
C--
C  Calculate chargino masses.
C--
      SQDET = AM2*AMU - AMW*AMW*SN2B
      B     = AM2*AM2 + AMU*AMU + 2*AMW*AMW
      D     = (B-2*SQDET)*(B+2*SQDET)
C>>>
      IF ( D.LT.0.D0 ) THEN
         RETURN
      ENDIF
C>>>
      SQD   = SQRT( D )
      AMX1  = SQRT( ( B - SQD )/2 )
      AMX2  = SQRT( ( B + SQD )/2 )
      AMXC(1) = AMX1
      AMXC(2) = AMX2
C--
C  Calculate mixing matrix.
C--
      PHIL  = ATAN( ( (AMX1-AM2)*(AMX1+AM2)-AMW*AMW*(1-CS2B) )
     .             /( SQ2*AMW*(AM2*CSB+AMU*SNB) ) )
      PHIL  = MOD(PHIL+PI,PI)
      PHIR  = ATAN( ( (AMX1-AM2)*(AMX1+AM2)-AMW*AMW*(1+CS2B) )
     .             /( SQ2*AMW*(AM2*SNB+AMU*CSB) ) )
      PHIR  = MOD(PHIR+PI,PI)
C--
      UL(1,1) = COS(PHIL)
      UL(1,2) = SIN(PHIL)
      UL(2,1) = -UL(1,2)
      UL(2,2) = UL(1,1)
C--
      UR(1,1) = COS(PHIR)
      UR(1,2) = SIN(PHIR)
      UR(2,1) = -UR(1,2)
      UR(2,2) = UR(1,1)
C--
      AMX1 = 0
      DO 100 J = 1, 2
         DO 10 I = 1, 2
            AMX1 = AMX1 + UR(1,I)*AMSMAT(I,J)*UL(1,J)
10       CONTINUE
100   CONTINUE
      IF ( AMX1.LT.0.D0 ) PHIR = PHIR + PI
C--
C  Decide EPSR.
C--
      EPSR = SIGN(1.D0,AM2*AMU-AMW*AMW*SN2B)
C--
C  That's it.
C--
      RETURN
      END
