C* (Update Record)
C*   91/09/02  K.Fujii       Reduced # bins for sig_0 tabulations.
C*
 
      SUBROUTINE SGAFBF(MODE,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                  E,SGEFF,ASEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE
      REAL*8       ALFS, ALF, S2W, AMZ, AMW, AMB, AMT, VTB2,
     .             E, SGEFF, ASEFF
      PARAMETER    ( NRSH = 1000 )
      INTEGER*4    NBN(3)
      REAL   *8    RSHRGN(0:3), DRSHRG(3),
     .             SGDAT(0:NRSH,3), ASDAT(0:NRSH,3)
      REAL   *4    XRS
      DATA         NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xPI    = ACOS(-1.D0)
         AME    = 0.511D-3
         ALF0   = 1./137.
         FACT   = 1 + 2*ALF0/xPI * ( xPI**2/6 - 1.D0/4 )
      ENDIF
C--
C  Set some variables.
C--
      IF ( MODE.LE.2 ) THEN
         IF ( AMT.LT.140.D0 ) THEN
            NBN   (1) = 5
            NBN   (2) = 160
            NBN   (3) = 10
            RSHRGN(0) = - 7.D0
            RSHRGN(1) = - 4.D0
            RSHRGN(2) =   0.D0
            RSHRGN(3) =   2.D0
            NBM       = 440
            DBM       = 1./NBM
         ELSE
            NBN   (1) = 5
            NBN   (2) = 70
            NBN   (3) = 10
            RSHRGN(0) = -12.D0
            RSHRGN(1) = - 7.D0
            RSHRGN(2) =   0.D0
            RSHRGN(3) =   5.D0
            NBM       = 100
            DBM       = 1./NBM
         ENDIF
      ENDIF
C--
C  Check if RS is in the acceptable range.
C--
      IF ( E.LT.RSHRGN(0) .OR. E.GT.RSHRGN(3) ) THEN
         PRINT *, ' SGAFBF: E = ', E, ' is out of range.'
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
         PRINT *, 'SGAFBF now tabulates sigma_0 for'
         PRINT *, '  ALFS   = ', ALFS
         PRINT *, '  AMT    = ', AMT
         PRINT *, '  VTB2   = ', VTB2
C>>>
         MODEP = MODE
         DO 100 IRGN = 1, 3
            MEH   = NBN(IRGN)
            EHMN  = RSHRGN(IRGN-1)
            EHMX  = RSHRGN(IRGN)
            DEH   = (EHMX-EHMN)/MEH
            DRSHRG(IRGN) = DEH
            DO 10 IEH = 0, MEH
               EH = EHMN + DEH*IEH
               CALL SGAFB0(MODEP,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                     EH, SIGTOT,DLTAFB)
               SGDAT(IEH,IRGN) = SIGTOT
               ASDAT(IEH,IRGN) = DLTAFB
C>>>
               PRINT *, EH, SIGTOT
C>>>
               MODEP = 3
10          CONTINUE
            MODEP = 2
100      CONTINUE
C>>>
         PRINT *, 'SGAFBF finished the tabulation of sigma_0 for'
C>>>
      ENDIF
C--
C  Gaussian beam width, beamstrahlung, and bremsstrahlng.
C--
      RS    = 2*AMT + E
      SGEFF = 0
      ASEFF = 0
      DO 20 IBM = 0, NBM
         IF ( IBM.EQ.0 .OR. IBM.EQ.NBM ) THEN
            WGT = 1.D0/3
         ELSE IF ( MOD(IBM,2).EQ.0 ) THEN
            WGT = 2.D0/3
         ELSE
            WGT = 4.D0/3
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
         MEH  = NBN(IRGN)
         EHMN = RSHRGN(IRGN-1)
         DEH  = DRSHRG(IRGN)
C--
         IEH  = (EH-EHMN)/DEH
         IF ( IEH.LT.MEH ) THEN
            F     = (EH-EHMN-IEH*DEH)/DEH
            SG    = SGDAT(IEH,IRGN)
     .              + (SGDAT(IEH+1,IRGN)-SGDAT(IEH,IRGN))*F
            AS    = ASDAT(IEH,IRGN)
     .              + (ASDAT(IEH+1,IRGN)-ASDAT(IEH,IRGN))*F
         ELSE
            SG    = SGDAT(MEH,IRGN)
            AS    = ASDAT(MEH,IRGN)
         ENDIF
         SGEFF = SGEFF + SG*WGT*DBM
         ASEFF = ASEFF + AS*SG*WGT*DBM
20    CONTINUE
C--
C  Vertex correction.
C--
      ASEFF = ASEFF/SGEFF
      SGEFF = FACT*SGEFF
C--
C  That's it.
C--
      RETURN
      END
