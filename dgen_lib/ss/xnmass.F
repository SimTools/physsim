 
      SUBROUTINE XNMASS(AM2,AMU,BT,SN2W,AMZ,AMXN,ON,ETA)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL   *8  AM2, AMU, BT, SN2W, AMZ, AMXN(4), ON(4,4)
      COMPLEX*16 ETA(4)
      REAL   *8  AMSMAT(10), AMX(4), EV(4,4), WORK(8)
      INTEGER*4  IP(4)
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.D0)
         SQ2   = SQRT(2.D0)
      ENDIF
C--
C  Calculate mass matrix.
C--
      CS2W = 1 - SN2W
      TN2W = SN2W/CS2W
      SNW  = SQRT(SN2W)
      CSW  = SQRT(CS2W)
      AM1  = (5.D0/3)*AM2*TN2W
      CSB  = COS(BT)
      SNB  = SIN(BT)
C--
      AMSMAT( 1) =  AM1
      AMSMAT( 2) =  0
      AMSMAT( 3) =  AM2
      AMSMAT( 4) = -AMZ*SNW*CSB
      AMSMAT( 5) =  AMZ*CSW*CSB
      AMSMAT( 6) =  0
      AMSMAT( 7) =  AMZ*SNW*SNB
      AMSMAT( 8) = -AMZ*CSW*SNB
      AMSMAT( 9) = -AMU
      AMSMAT(10) =  0
C--
C  Diagonalize mass matrix.
C--
      CALL DSEIG1(AMSMAT,4,AMX,EV,4,ME,WORK,IRT)
      IF ( IRT.NE.0 ) THEN
         PRINT *, ' >>>>> Error in XNMASS: IRT = ', IRT
      ENDIF
C--
C  Sort out mass eigen values.
C--
C>>>
C     IF ( NCALL.EQ.1 ) THEN
C        NCALL = 2
         CALL USORTD(4,1,1,AMX,IP)
C     ENDIF
C>>>
      DO 100 I = 1, 4
         AMXN(I) = ABS(AMX(IP(I)))
         IF ( AMX(IP(I)).GE.0.D0 ) THEN
            ETA(I) = 1
         ELSE
            ETA(I) = (0.D0,1.D0)
         ENDIF
         DO 10 J = 1, 4
            ON(I,J) = EV(J,IP(I))
10       CONTINUE
100   CONTINUE
C--
C  That's it.
C--
      RETURN
      END
