C* (Update Record)
C*   91/09/02  K.Fujii       Reduced # bins for sig_0 tabulations.
C*
 
      SUBROUTINE SGTTEF(MODE,AMT,ALFS,VTB2,AMH,BTH,RS,SGEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE
      REAL   *8    AMT, ALFS, VTB2, AMH, BTH, RS, SGEFF
      PARAMETER    ( NRSH = 100 )
      INTEGER*4    NBN(3)
      REAL   *8    RSHRGN(0:3), DRSHRG(3), SGDAT(0:NRSH,3)
      REAL   *4    XRS
      DATA NBN     /   5, 70,  10 /
      DATA RSHRGN  / -12.D0, -7.D0, 0.D0, 20.D0 /
      DATA       NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         xPI    = ACOS(-1.D0)
         S2W    = 0.230
         AME    = 0.511D-3
         AMB    = 5.
         AMZ    = 91.17
         AMW    = 80.
         GMZ    = 2.5
         ALF    = 1./128.
         ALF0   = 1./137.
         FACT   = 1 + 2*ALF0/xPI * ( xPI**2/6 - 1.D0/4 )
C>>>
         NBM    = 100
C        NBM    = 200
C>>>
         DBM    = 1./NBM
      ENDIF
C--
C  Set some variables.
C--
      VFF   = SQRT(VTB2)
C--
C  Check if RS is in the acceptable range.
C--
      E = RS - 2*AMT
      IF ( E.LT.RSHRGN(0) .OR. E.GT.RSHRGN(3) ) THEN
         PRINT *, ' SGTTEF: sqrt(s) = ', RS, ' is out of range.'
         PRINT *, ' AMT = ', AMT, ' EMN = ', RSHRGN(0),
     .            ' EMX = ', RSHRGN(3)
         PRINT *, '       : causes forced STOP.'
         STOP
      ENDIF
C--
C  Tabulate cross section w/o RC.
C--
      IF ( MODE.LE.2 ) THEN
C>>>
         PRINT *, 'SGTTEF now tabulates sigma_0 for'
         PRINT *, '  ALFS   = ', ALFS
         PRINT *, '  AMT    = ', AMT
         PRINT *, '  VTB2   = ', VFF**2
         PRINT *, '  AMH    = ', AMH
         PRINT *, '  BTH    = ', BTH
C>>>
         MODEP = MODE
         DO 100 IRGN = 1, 3
            MRSH  = NBN(IRGN)
            RSHMN = 2*AMT + RSHRGN(IRGN-1)
            RSHMX = 2*AMT + RSHRGN(IRGN)
            DRSH  = (RSHMX-RSHMN)/MRSH
            DRSHRG(IRGN) = DRSH
            DO 10 IRSH = 0, MRSH
               RSH = RSHMN + DRSH*IRSH
               CALL SGTTHR(MODEP,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,AMT,VFF,
     .                     AMH,BTH,   RSH,SG)
               SGDAT(IRSH,IRGN) = SG
               MODEP = 2
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
         EH  = RSH - 2*AMT
         IF ( EH.LE.RSHRGN(0) )                  GO TO 20
         IF ( EH.LE.RSHRGN(1) ) THEN
            IRGN = 1
         ELSE IF ( EH.LE.RSHRGN(2) ) THEN
            IRGN = 2
         ELSE
            IRGN = 3
         ENDIF
         MRSH  = NBN(IRGN)
         RSHMN = 2*AMT + RSHRGN(IRGN-1)
         DRSH  = DRSHRG(IRGN)
C--
         IRSH  = (RSH-RSHMN)/DRSH
         IF ( IRSH.LT.MRSH ) THEN
            F     = (RSH-RSHMN-IRSH*DRSH)/DRSH
            SG    = SGDAT(IRSH,IRGN)
     .              + (SGDAT(IRSH+1,IRGN)-SGDAT(IRSH,IRGN))*F
         ELSE
            SG    = SGDAT(MRSH,IRGN)
         ENDIF
         SGEFF = SGEFF + SG*WGT*DBM
20    CONTINUE
C--
C  Vertex correction.
C--
      SGEFF = FACT*SGEFF
C--
C  That's it.
C--
      RETURN
      END
