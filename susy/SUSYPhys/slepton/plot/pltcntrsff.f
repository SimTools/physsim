CINCLUDE SFMASS
C*
C*
C*
C*
C*
C*
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMU, GMSN, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMU, GMSN
C--
      INTEGER*4  I0(100)
      REAL   *8  AMSF2(7)
      PARAMETER ( MXxCT = 20 )
      REAL   *8  CONTUR(MXxCT), XYLIM(2,2), XY0(2), DXY0(2), CNT
      PARAMETER ( NP = 5000 )
      REAL   *8  XYDATA(3,0:NP)
      REAL   *8  FNCT
      EXTERNAL   FNCT
      DATA LOU   / 20 /
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C     ITYP = (1,2,3,4,5,6,7) = (SUL,SDL,SUR,SDR,SNL,SEL,SER)
C--
      PI   = ACOS(-1.)
      G2PB = 3.893796623D8
C--
      ITYP = 7
C      ITYP = 8
C      AMU  = 400
      AMU  = 600
C      TNBT =+2.
      TNBT =-2.
      BETA = MOD(ATAN(TNBT)+PI,PI)
C--
      S2W  = 0.23
      C2W  = 1 - S2W
      SNW  = SQRT(S2W)
      CSW  = SQRT(C2W)
      ALF  = 1/128.D0
      AMW  = 80.
      GMW  = 2.0
      AMZ  = 91.17
      GMZ  = 2.5
C--
C  Allocate output data set.
C--
C     CALL UALCPS(LOU,'TKSF.@.@TOP@.TDR','RENEW','FB',30,IRT)
C--
C  Initialize window parameters.
C--
      XYLIM(1,1) =    0.
      XYLIM(2,1) = 1200.
      XYLIM(1,2) =    0.
      XYLIM(2,2) = 1200.
C--
C  Initialize countour parameters.
C--
      NCT       = 4
      CONTUR(1) = 100
      CONTUR(2) = 250
      CONTUR(3) = 500
      CONTUR(4) = 750
C--
C      NCT       = 5
C      CONTUR(1) = 500
C      CONTUR(2) = 1000
C      CONTUR(3) = 1500
C      CONTUR(4) = 2000
C      CONTUR(5) = 2500
C--
C  Give a point within the innermost contour.
C     X0 = AMU
C     Y0 = AM2
C--
      AM00 =   10
      AM20 =   0
      DAM  =  XYLIM(2,2)/NP
      DXY0(1) = DAM
      DXY0(2) = DAM
C--
      DO 100 IP = 0, NP
         AM0 = AM00 + IP*DAM
         AM2 = AM20 + IP*DAM
         CALL SFMASS(AM0,AM2,AMU,BETA,ALF,S2W,AMZ,AMSF2)
         IF ( ITYP.LE.7 ) THEN
            AM = AMSF2(ITYP)
            AM = SIGN(SQRT(ABS(AM)),AM)
         ELSE IF ( ITYP.EQ.8 ) THEN
            AMS1 = SIGN(SQRT(ABS(AMSF2(1))),AMSF2(1))
            AMS2 = SIGN(SQRT(ABS(AMSF2(2))),AMSF2(2))
            AMS3 = SIGN(SQRT(ABS(AMSF2(3))),AMSF2(3))
            AMS4 = SIGN(SQRT(ABS(AMSF2(4))),AMSF2(4))
            IF ( AMS1.GT.0.D0 .AND. AMS2.GT.0.D0 .AND.
     .           AMS3.GT.0.D0 .AND. AMS4.GT.0.D0 ) THEN
               AM = ( AMS1 + AMS2 + AMS3 + AMS4 )/4
            ELSE
               AM = -999.
            ENDIF
         ENDIF
         XYDATA(1,IP) = AM0
         XYDATA(2,IP) = AM2
         XYDATA(3,IP) = AM
100   CONTINUE
C--
C  Draw contours.
C--
      DO 400 ICT = 1, NCT
         NFLP  = 0
         CNT   = CONTUR(ICT)
         DO 20 IP = 0, NP-1
            DELL = XYDATA(3,IP  ) - CNT
            DELR = XYDATA(3,IP+1) - CNT
            IF ( DELL*DELR.LT.0. .AND. ABS(DELL).LT.PI/4 .AND.
     .           ABS(DELR).LT.PI/4 ) THEN
               NFLP     = NFLP + 1
               IF ( DELL.LT.0. ) THEN
                  I0(NFLP) = IP
               ELSE
                  I0(NFLP) = IP + 1
               ENDIF
            ENDIF
20       CONTINUE
C--
C  Draw contours.
C--
         PRINT *, ' NFLP = ', NFLP
         IF ( NFLP.GT.0 ) THEN
            WRITE(LOU,*) '(                '
            WRITE(LOU,*) '( ******************************** '
            WRITE(LOU,*) '(    ICT  = ', ICT
            WRITE(LOU,*) '(    CNT  = ', CONTUR(ICT)
            WRITE(LOU,*) '(    NFLP = ', NFLP
            WRITE(LOU,*) '(                '
C--
C           WRITE(6  ,*) '(                '
C           WRITE(6  ,*) '( ******************************** '
C           WRITE(6  ,*) '(    ICT  = ', ICT
C           WRITE(6  ,*) '(    CNT  = ', CONTUR(ICT)
C           WRITE(6  ,*) '(    NFLP = ', NFLP
C           WRITE(6  ,*) '(                '
         ENDIF
         DO 300 IFLP = 1, NFLP
            IP = I0(IFLP)
C>>>
            WRITE(LOU,'('' ( IFLP,IP,X0,Y0,FL0,CNT = '',
     .                 I2,I7,2F10.3,2F8.5)')  IFLP, IP,
     .                 XYDATA(1,IP), XYDATA(2,IP), XYDATA(3,IP),
     .                 CONTUR(ICT)
C>>>
            XY0(1)  = XYDATA(1,IP)
            XY0(2)  = XYDATA(2,IP)
            CNT     = CONTUR(ICT)
            CALL UCONTR(FNCT,1,CONTUR(ICT),XYLIM,XY0,DXY0,LOU)
C>>>
C           WRITE(LOU,*) XY0(1), XY0(2)
C           WRITE(LOU,'('' SYMBOL 9O SIZE 0.5; PLOT'')')
C           WRITE(LOU,'('' SYMBOL 9O SIZE 0.1'')')
C>>>
300      CONTINUE
400   CONTINUE
C--
C  That's it.
C--
      STOP
      END
 
      DOUBLE PRECISION FUNCTION FNCT(AM0,AM2)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL   *8  AM0, AM2
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMU, GMSN, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMU, GMSN
C--
      REAL   *8  AMSF2(7)
      DATA NCALL /0/
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
C  Evaluate function value.
C--
      AM2P = AM2
      AM0P = AM0
      CALL SFMASS(AM0P,AM2P,AMU,BETA,ALF,S2W,AMZ,AMSF2)
C--
      IF ( ITYP.LE.7 ) THEN
         FNCT = SQRT(AMSF2(ITYP))
      ELSE IF ( ITYP.EQ.8 ) THEN
         FNCT = ( SQRT(AMSF2(1)) +
     .            SQRT(AMSF2(2)) +
     .            SQRT(AMSF2(3)) +
     .            SQRT(AMSF2(4)) )/4
      ELSE
         PRINT *, ' Invalid ITYP in FNCT. STOP.'
         STOP
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
