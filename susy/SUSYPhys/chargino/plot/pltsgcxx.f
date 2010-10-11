CINCLUDE SGCXXA
CINCLUDE DSGCXX
C*
C*  Total cross section for e+e- --> X+ X-.
C*    RS = sqrt(S)
C*
C*  5/16/91  K.Fujii       This version can handle X+1 X-2 case.
C*			   Tree level cross section.
C*
      IMPLICIT    REAL*8  ( A-H, O-Z )
      REAL   *8   GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2),
     .            AMSN, GMSN, AMZ, GMZ, AMX(2), GMX(2)
      REAL   *8   SG(0:3)
      REAL   *8   UL(2,2), UR(2,2)
C--
      CHARACTER*20 PATRN(0:3)
      DATA PATRN  / '0.2 0.0',
     .              '0.2 0.05 0.01 0.05',
     .              '0.2 0.05',
     .              '0.01 0.05 0.01' /
      CHARACTER*9 JOIN(0:3)
      DATA JOIN   / 'SOLID  ',
     .              'PATTERNED',
     .              'PATTERNED',
     .              'PATTERNED'/
C--
      DATA LOUNIT / 20 /
C
C========< Entry Point >================================================
C
C--
C  Total cross section.
C--
      WRITE(LOUNIT,'(''SET FONT DUPLEX'')')
      WRITE(LOUNIT,'(''NEW FRAME'')')
      WRITE(LOUNIT,'(''SET WINDOW X 2 10 Y 2 8'')')
      WRITE(LOUNIT,'(''SET LIMITS X 100. 1.E3 Y 0.001 10.'')')
      WRITE(LOUNIT,'(''SET SCALE Y LOGARITHMIC'')')
      WRITE(LOUNIT,'(''SET LABELS SIZE 3.0'')')
      WRITE(LOUNIT,'(''SET TITLE  SIZE 3.0'')')
      WRITE(LOUNIT,'(''TITLE   6. 7.5 SIZE 2.5'',
     .            '' '''' e2+3e2-3 --> C0615S262+3C0615S262-3'')')
      WRITE(LOUNIT,'(''CASE                   '',
     .            '' ''''  X X X X     GUUVVMVVX XGUUVVMVVX X'')')
      WRITE(LOUNIT,'(''TITLE 5.7 1.2 ''''sqrt(s)'')')
      WRITE(LOUNIT,'(''TITLE 0.2 3.7 ANGLE 90.'',
     .                            '' ''''S (pb)'')')
      WRITE(LOUNIT,'(''CASE          ''''G     '')')
      WRITE(LOUNIT,'(''SET ORDER X Y'')')
C--
C  Define constants.
C--
      x2PI      = 2*ACOS(-1.)
      x4PI      = 2*x2PI
C--
C  Standard Model parameters.
C--
      xSIN2W    = 0.23
      xCOS2W    = 1 - xSIN2W
      xSINW     = SQRT(xSIN2W)
      xCOSW     = SQRT(xCOS2W)
      xALF      = 1.D0/128D0
      AMZ       = 91.1
      GMZ       = 2.5
C--
C  Widths and masses.
C--
      IXM       = 1
      IXP       = 1
      GMX(1)    = 0.
      GMX(2)    = 0.
      AMSN      = 250.
      GMSN      = 0
C--
C  Chargino mass range.
C--
      NM   = 3
      AMMN = 100
      DM   = 50
C--
C  Polarization.
C--
      POL  = 0
C--
C  sqrt(s) range.
C--
      NRS  = 200
      RSMX = 1000
C--
C  Angular region.
C--
      CSMN = -1
      CSMX = -CSMN
      NCS  = 100
      DCS  = (CSMX-CSMN)/NCS
      PHI  = 0.
C--
C  Couplings of electron to Z and gamma.
C--
      C0        = SQRT(x4PI*xALF)
      CW        = C0/xSINW
      CZ        = CW/xCOSW
      QE        = -1
      T3LE      = -0.5D0
      T3RE      =  0
C--
      GAL(1)    = C0*QE
      GAL(2)    = C0*QE
      GZL(1)    = CZ*(T3LE-QE*xSIN2W)
      GZL(2)    = CZ*(T3RE-QE*xSIN2W)
C--
C  Loop over mixing angles.
C--
      NFI  = 2
      FIMN = 0
      FIMX = x2PI/4
      DFI  = (FIMX-FIMN)/NFI
      DO 10000 IFI = 0, NFI
C--
C  Coupling of chargino to Z and gamma.
C     U_i = ! cos(phi_i)  -sin(phi_i) !
C           ! sin(phi_i)   cos(phi_i) !
C  where i = (L,R). U_i transforms weak eigen states(wino,higgsino)
C  to mass eigen states(charginos):
C   ! X_- !   = !  cos(phi_i)  sin(phi_i) ! * ! W^- !
C   ! X_+ !_i   ! -sin(phi_i)  cos(phi_i) !   ! H^- !_i
C--
         FIL       = FIMN + DFI*IFI
         FIR       = FIL
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
         UR(2,1)   = -SNR
         UR(2,2)   =  CSR
C--
         QX        = -1
         T3LX      =   UL(IXM,1)*UL(IXP,1)*(-1)
     .               + UL(IXM,2)*UL(IXP,2)*(-1./2)
         T3RX      =   UR(IXM,1)*UR(IXP,1)*(-1)
     .               + UR(IXM,2)*UR(IXP,2)*(-1./2)
C--
         GAX(1)    = C0*QX*(UL(IXM,1)*UL(IXP,1)+UL(IXM,2)*UL(IXP,2))
         GAX(2)    = C0*QX*(UR(IXM,1)*UR(IXP,1)+UR(IXM,2)*UR(IXP,2))
         GZX(1)    = CZ*(T3LX-QX*xSIN2W)
         GZX(2)    = CZ*(T3RX-QX*xSIN2W)
         GSNX(1,1) = CW*UR(IXM,1)
         GSNX(2,1) = 0
         GSNX(1,2) = 0
         GSNX(2,2) = CW*UR(IXP,1)
         DO 1000 IM = 0, NM
            AM = AMMN + IM*DM
            RSMN = 2*AM + 5.E-2
            DRS  = (RSMX-RSMN)/NRS
            AMX(1) = AM
            AMX(2) = AM
            DO 100 IRS = 0, NRS
               RS   = RSMN + IRS*DRS
C>>>
               CALL SGCXXA(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .                    AMX, GMX ,RS, POL, SG)
               SGSM = SG(0)
C>>>
C                SGSM = 0
C                DO 10 ICS = 0, NCS
C                   IF ( ICS.EQ.0 .OR. ICS.EQ.NCS ) THEN
C                      WGT = 1./3
C                   ELSE IF ( MOD(ICS,2).EQ.0 ) THEN
C                      WGT = 2./3
C                   ELSE
C                      WGT = 4./3
C                   ENDIF
C                   CS = CSMN + DCS*ICS
C                   CALL DSGCXX(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
C      .                        AMX, GMX ,RS, POL, CS, PHI, SG)
C                   SGSM = SGSM + SG(0)*WGT
C 10             CONTINUE
C                PRINT *, ' RS, SGSM, SGA =', RS, SGSM*x2PI*DCS, SGSA
C>>>
               WRITE(LOUNIT,'(2E15.5)') RS, SGSM
100         CONTINUE
            WRITE(LOU,'(''SET PATTERN '',A20)') PATRN(IFI)
            WRITE(LOUNIT,'(''JOIN '',A9)') JOIN(IFI)
1000     CONTINUE
10000 CONTINUE
C--
      STOP
      END
