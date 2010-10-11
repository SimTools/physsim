C***********************************************************************
C*
C*===============================================================
C* Subrouine SGXXEF(MODE,GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
C*                                     AMX, GMX ,RS, POL, SGEFF)
C*========================================================-----==
C*
C* (Purpose)
C*   Calculates the effective total cross section with ISR and
C*   beam effects for chargino pair production.
C* (Inputs)
C*      MODE      : (I*4) : >2 skipps cross section tabulation.
C*      GAL(2)    : (R*4) : L-A-L chiral couplings.
C*      GAX(2)    : (R*4) : XC-A-XC chiral couplings.
C*      GZL(2)    : (R*4) : L-Z-L chiral couplings.
C*      GZX(2)    : (R*4) : XC-Z-XC chiral couplings.
C*      GSNX(2,2) : (R*4) : E-SNUE-XC couplings.
C*      AMSN      : (R*4) : sneutrino mass.
C*      GMSN      : (R*4) : sneutrino width (not used).
C*      AMZ       : (R*4) : Z mass.
C*      GMZ       : (R*4) : Z width.
C*      AMX(2)    : (R*4) : XC masses.
C*      GMX(2)    : (R*4) : XC widths (not used).
C*      RS        : (R*4) : sqrt(s).
C*      POL       : (R*4) : electron beam polarization.
C* (Output)
C*      SGEFF     : (R*4) : effective cross section.
C* (Relation)
C*   Invokes 
C*		RSHDIS, SGCXXA
C* (Update Record)
C*   91/05/16  K.Fujii       Original version.
C*   92/06/26  K.Fujii       Reduced # bins for sig_0 tabulations.
C*
C***********************************************************************
 
      SUBROUTINE SGXXEF(MODE,GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .                  AMX, GMX ,RS, POL, SGEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE
      REAL   *4    GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2), POL,
     .             AMSN, GMSN, AMZ, GMZ, AMX(2), GMX(2), RS, SGEFF
C--
      REAL   *8    GALD(2), GAXD(2), GZLD(2), GZXD(2), GSNXD(2,2),
     .             AMXD(2), GMXD(2), SG(0:3)
C--
      REAL   *4    XRS
C--
      PARAMETER    ( NRSH = 200 )
      INTEGER*4    NBN(3)
      DATA NBN     /  50, 20, 100 /
      REAL   *8    RSHRGN(0:3), DRSHRG(3), SGDAT(0:NRSH,3)
      DATA RSHRGN  /   0.D0,  5.D0, 25.D0, 1025.D0 /
      DATA       NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NBM    = 100
         DBM    = 1./NBM
      ENDIF
C--
C  Check if RS is in the acceptable range.
C--
      E = RS - AMX(1) - AMX(2)
      IF ( E.LT.RSHRGN(0) .OR. E.GT.RSHRGN(3) ) THEN
         PRINT *, ' SGTTEF: sqrt(s) = ', RS, ' is out of range.'
         PRINT *, ' AMX = ', AMX, ' EMN = ', RSHRGN(0),
     .            ' EMX = ', RSHRGN(3)
         PRINT *, '       : causes forced STOP.'
         STOP
      ENDIF
C--
C  Tabulate cross section w/o RC.
C--
      IF ( MODE.LE.2 ) THEN
         PRINT *, 'SGXXEF now tabulates sigma_0.'
C--
         GALD(1)    = GAL(1)
         GALD(2)    = GAL(2)
         GAXD(1)    = GAX(1)
         GAXD(2)    = GAX(2)
         GZLD(1)    = GZL(1)
         GZLD(2)    = GZL(2)
         GZXD(1)    = GZX(1)
         GZXD(2)    = GZX(2)
         GSNXD(1,1) = GSNX(1,1)
         GSNXD(2,1) = GSNX(2,1)
         GSNXD(1,2) = GSNX(1,2)
         GSNXD(2,2) = GSNX(2,2)
         AMSND      = AMSN
         GMSND      = GMSN
         AMZD       = AMZ
         GMZD       = GMZ
         AMXD(1)    = AMX(1)
         AMXD(2)    = AMX(2)
         GMXD(1)    = GMX(1)
         GMXD(2)    = GMX(2)
         AMX12      = AMX(1) + AMX(2)
         POLD       = POL
C>>>
C        PRINT *, ' GALD       = ',  GALD
C        PRINT *, ' GAXD       = ',  GAXD
C        PRINT *, ' GAXD       = ',  GAXD
C        PRINT *, ' GZLD       = ',  GZLD
C        PRINT *, ' GZLD       = ',  GZLD
C        PRINT *, ' GZXD       = ',  GZXD
C        PRINT *, ' GZXD       = ',  GZXD
C        PRINT *, ' GSNXD(*,1) = ',  GSNXD(1,1), GSNXD(2,1)
C        PRINT *, ' GSNXD(*,2) = ',  GSNXD(1,2), GSNXD(2,2)
C        PRINT *, ' AMSND      = ',  AMSND
C        PRINT *, ' GMSND      = ',  GMSND
C        PRINT *, ' AMZD       = ',  AMZD
C        PRINT *, ' GMZD       = ',  GMZD
C        PRINT *, ' AMXD       = ',  AMXD
C        PRINT *, ' GMXD       = ',  GMXD
C        PRINT *, ' GMXD       = ',  AMX12
C>>>
C--
         DO 100 IRGN = 1, 3
            MRSH  = NBN(IRGN)
            RSHMN = AMX12 + RSHRGN(IRGN-1)
            RSHMX = AMX12 + RSHRGN(IRGN)
            DRSH  = (RSHMX-RSHMN)/MRSH
            DRSHRG(IRGN) = DRSH
            DO 10 IRSH = 0, MRSH
               RSH = RSHMN + DRSH*IRSH
               CALL SGCXXA(GALD,GAXD,GZLD,GZXD,GSNXD, AMSND,GMSND,
     .                     AMZD,GMZD,AMXD,GMXD, RSH, POLD, SG)
               SGDAT(IRSH,IRGN) = SG(0)
10          CONTINUE
100      CONTINUE
      ENDIF
C--
C  Gaussian beam width, beamstrahlung, and bremsstrahlng.
C--
      SGEFF = 0
      DO 20 IBM = 0, NBM
         IF ( IBM.EQ.0 .OR. IBM.EQ.NBM ) THEN
            WGT = 1./3
         ELSE IF ( MOD(IBM,2).EQ.0 ) THEN
            WGT = 2./3
         ELSE
            WGT = 4./3
         ENDIF
         CALL RSHDIS(REAL(DBM*IBM),1,XRS)
         RSH = RS*XRS
C--
         EH  = RSH - AMX12
         IF ( EH.LE.RSHRGN(0) )                  GO TO 20
         IF ( EH.LE.RSHRGN(1) ) THEN
            IRGN = 1
         ELSE IF ( EH.LE.RSHRGN(2) ) THEN
            IRGN = 2
         ELSE
            IRGN = 3
         ENDIF
         MRSH  = NBN(IRGN)
         RSHMN = AMX12 + RSHRGN(IRGN-1)
         DRSH  = DRSHRG(IRGN)
C--
         IRSH  = (RSH-RSHMN)/DRSH
         IF ( IRSH.LT.MRSH ) THEN
            F     = (RSH-RSHMN-IRSH*DRSH)/DRSH
            SG0   = SGDAT(IRSH,IRGN)
     .              + (SGDAT(IRSH+1,IRGN)-SGDAT(IRSH,IRGN))*F
         ELSE
            SG0   = SGDAT(MRSH,IRGN)
         ENDIF
         SGEFF = SGEFF + SG0*WGT*DBM
20    CONTINUE
C--
C  That's it.
C--
      RETURN
      END
