CINCLUDE XCMASS
CINCLUDE SGCXXA
C*
C*  Contour plots for charginos.
C*
C*  5/16/91  K.Fujii       Original version.
C*
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSN, GMSN, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSN, GMSN
C--
      INTEGER*4  I0(100)
      REAL   *8  AMXC(2)
C>>>
C      PARAMETER ( NCT = 25 )
      PARAMETER ( NCT = 4 )
C>>>      
      REAL   *8  CONTUR(NCT), XYLIM(2,2), XY0(2), DXY0(2), CNT
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
C--
      PI   = ACOS(-1.)
      G2PB = 3.893796623D8
C--
C  Job parameters.
C     ITYP = 1 : mixing angle phi_L
C	   = 2 : mixing angle phi_R
C          = 3 : lighter chargino mass chi1
C          = 4 : heavier chargino mass chi2
C          = 5 : sig_tot for chi1-chi1
C          = 6 : sig_tot for chi1-chi2
C          = 7 : sig_tot for chi2-chi2
C     RS   = sqrts in GeV
C     TNBT = tan(beta)
C     AMSN = electron-sneutrino mass
C--
      ITYP = 3
      RS   = 500.
C      TNBT =+2.
      TNBT =-2.
      BETA = MOD(ATAN(TNBT)+PI,PI)
      AMSN = 500.
C--
C  Standard Model parameters.
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
C--
C  Allocate output data set.
C--
C     CALL UALCPS(LOU,'TKSF.@.@TOP@.TDR','RENEW','FB',30,IRT)
C--
C  Initialize window parameters.
C--
      IF ( ITYP.LE.4 ) THEN
         XYLIM(1,1) =    0.
         XYLIM(2,1) = 1000.
         XYLIM(1,2) =    0.
         XYLIM(2,2) = 1000.
      ELSE IF ( ITYP.GE.5 .AND. ITYP.LE.7 ) THEN
         XYLIM(1,1) =    0.
         XYLIM(2,1) = 1000.
         XYLIM(1,2) =    0.
         XYLIM(2,2) = 1000.
      ENDIF
C--
C  Initialize countour parameters.
C--
      LOU = 20
      IF ( ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
         DFI = 2*PI/NCT
      ELSE IF ( ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
         DFI = XYLIM(2,1)/NCT
      ELSE IF ( ITYP.EQ.5 .OR. ITYP.EQ.6 .OR. ITYP.EQ.7 ) THEN
         DFI = 2.5D0/NCT
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
      IF ( ITYP.EQ.1 ) THEN
C--For phi_L
         AM20 =   0.
         AMU0 =  XYLIM(2,2)/2
         DAM  =  AMU0/NP
         DXY0(1) = MAX(DAM,200.)
         DXY0(2) = MAX(DAM,200.)
      ELSE IF ( ITYP.EQ.2 ) THEN
C--For phi_R
         AM20 =   0.
         AMU0 =  XYLIM(2,2)
         DAM  =  AMU0/NP
C        DXY0(1) = DAM/10
C        DXY0(2) = DAM/10
         DXY0(1) = DAM
         DXY0(2) = DAM
      ELSE IF ( ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
C--For masses.
         AMU0 =   1
         AM20 =   0
         DAM  =  XYLIM(2,2)/NP
         DXY0(1) = DAM
         DXY0(2) = DAM
      ELSE IF ( ITYP.GE.5 .AND. ITYP.LE.7 ) THEN
C--For cross sections.
         AMU0 =    0
         AM20 =    5
C        AMU0 =    0
C        AM20 =    900
         DAM  =  XYLIM(2,2)/NP
         DAMU =  DAM*1
         DAM2 =  DAM*0
         DXY0(1) = MAX(DAM,100.)
         DXY0(2) = MAX(DAM,100.)
      ENDIF
C--
      DO 100 IP = 0, NP
         IF ( ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            AMU = AMU0 - IP*DAM
            AM2 = AM20 + IP*DAM
         ELSE IF ( ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            AMU = AMU0 + IP*DAM
            AM2 = AM20 + IP*DAM
         ELSE IF ( ITYP.GE.5 .AND. ITYP.LE.7 ) THEN
            AMU = AMU0 + IP*DAMU
            AM2 = AM20 + IP*DAM2
         ENDIF
C>>>
         CALL XCMASS(AM2,AMU,BETA,AMW,AMXC,PHIL,PHIR,EPSR)
         IF ( ITYP.EQ.2 ) THEN
            PHIL = PHIR
         ELSE IF ( ITYP.EQ.3 ) THEN
            PHIL = AMXC(1)
         ELSE IF ( ITYP.EQ.4 ) THEN
            PHIL = AMXC(2)
         ELSE IF ( ITYP.GE.5 .AND. ITYP.LE.7 ) THEN
            CALL SGTOT(AM2,AMU,RRAT)
C>>>
         IF ( RRAT.GT.0.D0 ) THEN
C           PRINT 66,  AMU,AM2, PHIL,PHIR, EPSR, AMXC(1),AMXC(2), RRAT
66          FORMAT(2F7.2, 2F8.4, F6.1, 2F8.2, F15.5)
         ENDIF
C>>>
            PHIL = RRAT
         ENDIF
C>>>
         XYDATA(1,IP) = AMU
         XYDATA(2,IP) = AM2
         XYDATA(3,IP) = PHIL
C>>>
C        PRINT *, ' AMU, AM2, PHIL = ', AMU, AM2, PHIL
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
 
      DOUBLE PRECISION FUNCTION FNCT(AMU,AM2)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL   *8  AMU, AM2
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSN, GMSN, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSN, GMSN
C--
      REAL   *8  AMXC(2)
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
      IF ( ITYP.LE.4 )
     .     CALL XCMASS(AM2P,AMUP,BETA,AMW,AMXC,PHIL,PHIR,EPSR)
C--
      IF ( ITYP.EQ.1 ) THEN
         FNCT = PHIL
      ELSE IF ( ITYP.EQ.2 ) THEN
         FNCT = PHIR
      ELSE IF ( ITYP.EQ.3 ) THEN
         FNCT = AMXC(1)
      ELSE IF ( ITYP.EQ.4 ) THEN
         FNCT = AMXC(2)
      ELSE IF ( ITYP.GE.5 .AND. ITYP.LE.7 ) THEN
         CALL SGTOT(AM2P,AMUP,RRAT)
         FNCT = RRAT
      ELSE
         PRINT *, ' Invalid ITYP in FNCT. STOP.'
         STOP
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
 
      SUBROUTINE SGTOT(AM2,AMU,RRAT)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  AM2, AMU, RRAT
C--
      COMMON     AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSN, GMSN, ITYP
      INTEGER*4  ITYP
      REAL   *8  AMW, GMW, AMZ, GMZ, BETA, RS, S2W, C2W,
     .           SNW, CSW, ALF, G2PB, AMSN, GMSN
C--
      REAL   *8   GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2),
     .            AMXC(2), AMX(2), GMX(2)
      REAL   *8   SG(0:3)
      REAL   *8   UL(2,2), UR(2,2)
      DATA NCALL  /0/
C
C========< Entry Point >================================================
C
C--
C  Define constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         x2PI = 2*ACOS(-1.)
         x4PI = 2*x2PI
         SGPT = x4PI*ALF*ALF*G2PB/3
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
      ENDIF
C--
C  Widths and masses.
C--
      CALL XCMASS(AM2,AMU,BETA,AMW,AMXC,FIL,FIR,EPSR)
C--
      AMX(1)    = AMXC(IXM)
      AMX(2)    = AMXC(IXP)
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
      GAL(1)    = C0*QE
      GAL(2)    = C0*QE
      GZL(1)    = CZ*(T3LE-QE*S2W)
      GZL(2)    = CZ*(T3RE-QE*S2W)
C--
C  Coupling of chargino to Z and gamma.
C     U_i = ! cos(phi_i)  -sin(phi_i) !
C           ! sin(phi_i)   cos(phi_i) !
C  where i = (L,R). U_i transforms weak eigen states(wino,higgsino)
C  to mass eigen states(charginos):
C   ! X_- !   = !  cos(phi_i)  sin(phi_i) ! * ! W^- !
C   ! X_+ !_i   ! -sin(phi_i)  cos(phi_i) !   ! H^- !_i
C--
      CSL       = COS(FIL)
      SNL       = SIN(FIL)
      CSR       = COS(FIR)
      SNR       = SIN(FIR)
C--
      UL(1,1)   =  CSL
      UL(1,2)   =  SNL
      UL(2,1)   = -SNL
      UL(2,2)   =  CSL
C--
      UR(1,1)   =  CSR
      UR(1,2)   =  SNR
      UR(2,1)   = -SNR*EPSR
      UR(2,2)   =  CSR*EPSR
C--
      QX        =  -1
      T3LX      =   UL(IXM,1)*UL(IXP,1)*(-1)
     .            + UL(IXM,2)*UL(IXP,2)*(-1./2)
      T3RX      =   UR(IXM,1)*UR(IXP,1)*(-1)
     .            + UR(IXM,2)*UR(IXP,2)*(-1./2)
C--
      GAX(1)    = C0*QX*(UL(IXM,1)*UL(IXP,1)+UL(IXM,2)*UL(IXP,2))
      GAX(2)    = C0*QX*(UR(IXM,1)*UR(IXP,1)+UR(IXM,2)*UR(IXP,2))
      GZX(1)    = CZ*(T3LX-QX*S2W)
      GZX(2)    = CZ*(T3RX-QX*S2W)
      GSNX(1,1) = CW*UR(IXM,1)
      GSNX(2,1) = 0
      GSNX(1,2) = 0
      GSNX(2,2) = CW*UR(IXP,1)
C--
C  Calculate total cross section.
C--
      POL  = 0
      CALL SGCXXA(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .           AMX, GMX ,RS, POL, SG)
      RRAT = SG(0)*RS*RS/SGPT
C--
C  That's it.
C--
      RETURN
      END
