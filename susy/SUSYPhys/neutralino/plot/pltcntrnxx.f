CINCLUDE XNMASS
CINCLUDE SFMASS
CINCLUDE SGNXXA
CINCLUDE USORTD
C*
C*
C*
C*
C*
C*
      IMPLICIT   REAL*8 ( A-H, O-Z  )
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSE, GMSE, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSE(2), GMSE(2)
C--
      INTEGER*4  I0(100)
      REAL   *8  AMXN(4), ON(4,4)
      COMPLEX*16 ETA(4)
      PARAMETER ( NCT = 10 )
      REAL   *8  CONTUR(NCT), XYLIM(2,2), XY0(2), DXY0(2), CNT
      PARAMETER ( NP = 10000 )
      REAL   *8  XYDATA(3,0:NP)
      REAL   *8  FNCT
      EXTERNAL   FNCT
      DATA LOU   / 20 /
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      PI   = ACOS(-1.)
      G2PB = 3.893796623D8
C--
      ITYP = 1
      RS   = 1000.
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
      AMZ  = 91.1
      GMZ  = 2.5
      AMSE(2) = 500
C--
C  Allocate output data set.
C--
C      CALL UALCPS(LOU,'TKSF.@.@TOP@.TDR','RENEW','FB',30,IRT)
C--
C  Initialize window parameters.
C--
      IF ( ITYP.LE.4 ) THEN
         XYLIM(1,1) =    0.
         XYLIM(2,1) = 1000.
         XYLIM(1,2) =    0.
         XYLIM(2,2) = 1000.
      ELSE
         XYLIM(1,1) =    0.
         XYLIM(2,1) = 1000.
         XYLIM(1,2) =    0.
         XYLIM(2,2) = 1000.
      ENDIF
C--
C  Initialize countour parameters.
C--
      LOU = 20
      IF ( ITYP.LE.4 ) THEN
         DFI = XYLIM(2,1)/NCT
      ELSE
         DFI = 1.0D0/NCT
      ENDIF
C--
      DO 10 IFI = 1, NCT
         CONTUR(IFI) = IFI*DFI
10    CONTINUE
C--
C  Give a point within the innermost contour.
C     X0 = AMU
C     Y0 = AM2
C--
      IF ( ITYP.LE.4 ) THEN
C--For masses.
         AMU0 =   1
         AM20 =   0
         DAM  =  XYLIM(2,2)/NP
         DXY0(1) = DAM
         DXY0(2) = DAM
      ELSE
C--For cross sections.
         DAM  =  XYLIM(2,2)/NP
C        AMU0 =    0
C        AM20 =   1000
C        DAMU =  DAM*1
C        DAM2 = -DAM*0
C--
         AMU0 =   1000
         AM20 =    0
         DAMU = -DAM*1
         DAM2 =  DAM*1
         DXY0(1) = MAX(DAM,100.)
         DXY0(2) = MAX(DAM,100.)
      ENDIF
C--
      DO 100 IP = 0, NP
         IF ( ITYP.LE.4 ) THEN
            AMU = AMU0 + IP*DAM
            AM2 = AM20 + IP*DAM
         ELSE
            AMU = AMU0 + IP*DAMU
            AM2 = AM20 + IP*DAM2
         ENDIF
         CALL XNMASS(AM2,AMU,BETA,S2W,AMZ,AMXN,ON,ETA)
         IF ( ITYP.LE.4 ) THEN
            PHIL = AMXN(ITYP)
         ELSE
            CALL SGTOT(AM2,AMU,RRAT)
            PHIL = RRAT
         ENDIF
         XYDATA(1,IP) = AMU
         XYDATA(2,IP) = AM2
         XYDATA(3,IP) = PHIL
C>>>
C        PRINT *, ' AMU, AM2, SG = ', AMU, AM2, PHIL
C>>>
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
            WRITE(LOU,'('' ( IFLP,IP,X0,Y0,FL0,CNT = '',
     .                 I2,I7,2F10.3,2F8.5)')  IFLP, IP,
     .                 XYDATA(1,IP), XYDATA(2,IP), XYDATA(3,IP),
     .                 CONTUR(ICT)
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
 
      DOUBLE PRECISION FUNCTION FNCT(AMU,AM2)
 
      IMPLICIT   REAL*8 ( A-H, O-Z  )
      REAL   *8  AMU, AM2
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSE, GMSE, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSE(2), GMSE(2)
C--
      REAL   *8  AMXN(4), ON(4,4)
      COMPLEX*16 ETA(4)
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
      AMUP = AMU
      IF ( ITYP.LE.4 ) CALL XNMASS(AM2P,AMUP,BETA,S2W,AMZ,AMXN,ON,ETA)
C--
      IF ( ITYP.LE.4 ) THEN
         FNCT = AMXN(ITYP)
      ELSE
         CALL SGTOT(AM2P,AMUP,RRAT)
         FNCT = RRAT
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
 
C* 93/03/19   Factor 2 is corrected for ZXX vertices (TTFT).
C*
 
      SUBROUTINE SGTOT(AM2,AMU,RRAT)
 
      IMPLICIT   REAL*8  ( A-H, O-Y  )
      IMPLICIT   COMPLEX*16 ( Z )
      REAL   *8  AM2, AMU, RRAT
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSE, GMSE, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSE(2), GMSE(2)
C--
      REAL   *8   GZL(2), AMXN(4), AMX(2), GMX(2), AMSF(7)
      COMPLEX*16  GZX(2), GSEX(2,2)
      REAL   *8   SG(0:6)
      REAL   *8   ON(4,4)
      COMPLEX*16  ETA(4)
      DATA NCALL  /0/
C
C========< Entry Point >================================================
C
C--
C  Define constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         x2PI  = 2*ACOS(-1.)
         x4PI  = 2*x2PI
         SQ2   = SQRT(2.D0)
         SGPT  = x4PI*ALF*ALF*G2PB/3
         AM0   = 0
      ENDIF
C--
C  Decide reaction type.
C--
      IF ( ITYP.EQ.5 ) THEN
         IXM       = 1
         IXP       = 1
      ELSE IF ( ITYP.EQ.6 ) THEN
         IXM       = 1
         IXP       = 2
      ELSE IF ( ITYP.EQ.7 ) THEN
         IXM       = 2
         IXP       = 2
      ELSE IF ( ITYP.EQ.8 ) THEN
         IXM       = 1
         IXP       = 3
      ELSE IF ( ITYP.EQ.9 ) THEN
         IXM       = 2
         IXP       = 3
      ELSE IF ( ITYP.EQ.10 ) THEN
         IXM       = 3
         IXP       = 3
      ELSE IF ( ITYP.EQ.11 ) THEN
         IXM       = 1
         IXP       = 4
      ELSE IF ( ITYP.EQ.12 ) THEN
         IXM       = 2
         IXP       = 4
      ELSE IF ( ITYP.EQ.13 ) THEN
         IXM       = 3
         IXP       = 4
      ELSE IF ( ITYP.EQ.14 ) THEN
         IXM       = 4
         IXP       = 4
      ENDIF
C--
C  Widths and masses.
C--
      CALL XNMASS(AM2,AMU,BETA,S2W,AMZ,AMXN,ON,ETA)
      CALL SFMASS(AM0,AM2,AMU,BETA,ALF,S2W,AMZ,AMSF)
      AMSE(1) = SQRT( AMSE(2)**2 - AMSF(7) + AMSF(5) )
C--
      AMX(1)    = AMXN(IXM)
      AMX(2)    = AMXN(IXP)
      GMX(1)    = 0.
      GMX(2)    = 0.
C--
      QVAL = RS - AMX(1) - AMX(2)
      IF ( QVAL.LE.0.D0 ) THEN
         RRAT = 0
         RETURN
      ENDIF
C--
C  Couplings of electron to Z and gamma.
C--
      C0        = SQRT(x4PI*ALF)
      CW        = C0/SNW
      CZ        = CW/CSW
      QE        = -1
      T3LE      = -0.5D0
      T3RE      =  0
C--
      GZL(1)    = CZ*(T3LE-QE*S2W)
      GZL(2)    = CZ*(T3RE-QE*S2W)
C--
C  Coupling of neutralino to Z.
C--
      ZV        =  CZ/2*(0.D0,1.D0)*IMAG(ETA(IXM)*CONJG(ETA(IXP)))
     .             *( ON(IXM,3)*ON(IXP,3)-ON(IXM,4)*ON(IXP,4) )
      ZA        = -CZ/2            *DBLE(ETA(IXM)*CONJG(ETA(IXP)))
     .             *( ON(IXM,3)*ON(IXP,3)-ON(IXM,4)*ON(IXP,4) )
      IF ( IXM.EQ.IXP ) ZA = ZA/2
C--
      GZX(1)    = ZV - ZA
      GZX(2)    = ZV + ZA
      GSEX(1,1) = SQ2*CZ*ETA(IXM)
     .             *(T3LE*CSW*ON(IXM,2)+(QE-T3LE)*SNW*ON(IXM,1))
      GSEX(2,1) = SQ2*CZ*ETA(IXP)
     .             *(T3LE*CSW*ON(IXP,2)+(QE-T3LE)*SNW*ON(IXP,1))
      GSEX(1,2) = -SQ2*CZ*CONJG(ETA(IXM))*QE*SNW*ON(IXM,1)
      GSEX(2,2) = -SQ2*CZ*CONJG(ETA(IXP))*QE*SNW*ON(IXP,1)
C--
C  Convert GSEX convention to match SUSY manual.
C--
      GSEX(1,1) = CONJG(GSEX(1,1))
      GSEX(2,1) = CONJG(GSEX(2,1))
      GSEX(1,2) = CONJG(GSEX(1,2))
      GSEX(2,2) = CONJG(GSEX(2,2))
C--
C  Calculate total cross section.
C--
C     IF ( NCALL.EQ.1 ) THEN
C        NCALL = 2
C        PRINT *, ' AMX (  1) = ', AMX(1)
C        PRINT *, ' AMX (  2) = ', AMX(2)
C        PRINT *, ' GZX (  1) = ', GZX(1)
C        PRINT *, ' GZX (  2) = ', GZX(2)
C        PRINT *, ' GSEX(1,1) = ', GSEX(1,1)
C        PRINT *, ' GSEX(2,1) = ', GSEX(2,1)
C        PRINT *, ' GSEX(1,2) = ', GSEX(1,2)
C        PRINT *, ' GSEX(2,2) = ', GSEX(2,2)
C     ENDIF
      POL = 0
      CALL SGNXXA(GZL,GZX,GSEX,AMSE,GMSE,AMZ,GMZ,AMX,GMX,RS,POL,SG)
      RRAT = SG(0)*RS*RS/SGPT
      IF ( IXM.EQ.IXP ) RRAT = RRAT/2
C--
C  That's it.
C--
      RETURN
      END
