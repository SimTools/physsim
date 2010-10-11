CX*INCLUDE RSHDIS
CINCLUDE RSHDISN
CINCLUDE UDSRCH
C--
CINCLUDE COUPLSUB
CINCLUDE SFMASS
CINCLUDE INOMIX
CINCLUDE SWPELM
CINCLUDE INIPRM
CINCLUDE SGSESE
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  AM(2), AMX(4), CGE(2), CZE(2), CGSE(2), CZSE(2)
      COMPLEX*16 CX(2,2,4)
      DATA LOU   /20/
C
C==================< Entry Point >======================================
C
C--
C  Set experimental parameters.
C--
C      POLE = +1
      POLE = +0.9
      RSMX = 350
C--
C  Set MSSM parameters.
C--
      IFLG  = 1
C     AM0   = 150.D0
C     AMU   = 400.D0
C     AM2   = 283.D0
C     TNB   = 2.D0
C--
      AM0   =  70.D0
      AMU   = 400.D0
      AM2   = 250.D0
C      TNB   = 2.D0
      TNB   = 3.D0
      ITYP  = 2
C--
      x2PI  = 2*ACOS(-1.D0)
C--
      CALL INIPRM(4,IFLG,AM0,AMU,AM2,TNB,
     .            AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE)
C--
      WRITE(LOU,'(''( AM0  = '',F8.2)') AM0
      WRITE(LOU,'(''( AMU  = '',F8.2)') AMU
      WRITE(LOU,'(''( AM2  = '',F8.2)') AM2
      WRITE(LOU,'(''( TNB  = '',F8.2)') TNB
      WRITE(LOU,'(''( AMEL = '',F8.2)') AM(1)
      WRITE(LOU,'(''( AMER = '',F8.2)') AM(2)
      WRITE(LOU,'(''( AMX1 = '',F8.2)') AMX(1)
      WRITE(LOU,'(''( AMX2 = '',F8.2)') AMX(2)
      WRITE(LOU,'(''( AMX3 = '',F8.2)') AMX(3)
      WRITE(LOU,'(''( AMX4 = '',F8.2)') AMX(4)
C--
C  Loop over final state chirality combinations.
C--
      AMEL = AM(1)
      AMER = AM(2)
C--
C  Set selectron masses.
C--
      IF ( ITYP.EQ.1 ) THEN
         AM(1) = AMEL
         AM(2) = AMEL
      ELSE IF ( ITYP.EQ.2 ) THEN
         AM(1) = AMER
         AM(2) = AMER
      ELSE IF ( ITYP.EQ.3 ) THEN
         AM(1) = AMER
         AM(2) = AMEL
      ELSE IF ( ITYP.EQ.4 ) THEN
         AM(1) = AMEL
         AM(2) = AMER
      ENDIF
C--
C  Prepare topdraw data.
C--
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET WINDOW X 2.5 12.5 Y 2.2 9.6'')')
      WRITE(LOU,'(''SET LIMITS X 270 500 Y 5.E-4 5.E+0'')')
      WRITE(LOU,'(''SET SCALE Y LOG                   '')')
      WRITE(LOU,'(''SET LABELS SIZE 3      '')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''TITLE 6.4 1.2 SIZE 5 ''''2s0O (GeV)'')')
      WRITE(LOU,'(''CASE                 ''''M UD      '')')
      WRITE(LOU,'(''TITLE 0.4 4.7 SIZE 5 ANGLE 90 ''''S0tot1 (pb) '')')
      WRITE(LOU,'(''CASE                          ''''GX   X      '')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET ORDER X DUMMY Y '')')
C--
C  Loop over sqrt(s).
C--
      NRS  = 100
      RSMN = AM(1) + AM(2) + 1.D-3
      DRS  = (RSMX-RSMN)/NRS
      CSMX = 1.D0
C--
      MODE = 0
      DO 100 IRS = 0, NRS
         RS = RSMN + IRS*DRS
         CALL SGSESE(ITYP,AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE,
     .               RS,POLE, CSMX, SG0)
         CALL SGSESB(MODE,ITYP,AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE,
     .               RS,POLE, CSMX, SG)
C        PRINT *, ' RS = ', RS, ' GeV  sig  = ', SG , ' pb',
C    .                               ' sig0 = ', SG0, ' pb'
C--
         WRITE(LOU,'(3E15.5)') RS, SG0, SG
         MODE = 3
100   CONTINUE
      WRITE(LOU,'(''JOIN SOLID'')')
C--
C  That's it.
C--
      STOP
      END
C*
C* (Update Record)
C*   92/06/26  K.Fujii       Reduced # bins for sig_0 tabulations.
C*
 
      SUBROUTINE SGSESB(MODE,ITYP,AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE,
     .                  RS,POLE,CSMX, SGEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE, ITYP
      REAL   *8    RS, POLE, CSMX, AMZ, AM(2), AMX(4),
     .             CGE(2), CZE(2), CGSE(2), CZSE(2), SGEFF
      COMPLEX*16   CX(2,2,4)
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
      E = RS - AM(1) - AM(2)
      IF ( E.LT.RSHRGN(0) .OR. E.GT.RSHRGN(3) ) THEN
         PRINT *, ' SGTTEF: sqrt(s) = ', RS, ' is out of range.'
         PRINT *, ' AM = ', AM, ' EMN = ', RSHRGN(0),
     .            ' EMX = ', RSHRGN(3)
         PRINT *, '       : causes forced STOP.'
         STOP
      ENDIF
C--
C  Tabulate cross section w/o RC.
C--
      IF ( MODE.LE.2 ) THEN
         PRINT *, 'SGSESB now tabulates sigma_0.'
         AM12       = AM(1) + AM(2)
C--
         DO 100 IRGN = 1, 3
            MRSH  = NBN(IRGN)
            RSHMN = AM12 + RSHRGN(IRGN-1)
            RSHMX = AM12 + RSHRGN(IRGN)
            DRSH  = (RSHMX-RSHMN)/MRSH
            DRSHRG(IRGN) = DRSH
            DO 10 IRSH = 0, MRSH
               RSH = RSHMN + DRSH*IRSH
               CALL SGSESE(ITYP,AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE,
     .               RSH,POLE, CSMX, SG0)
               SGDAT(IRSH,IRGN) = SG0
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
         EH  = RSH - AM12
         IF ( EH.LE.RSHRGN(0) )                  GO TO 20
         IF ( EH.LE.RSHRGN(1) ) THEN
            IRGN = 1
         ELSE IF ( EH.LE.RSHRGN(2) ) THEN
            IRGN = 2
         ELSE
            IRGN = 3
         ENDIF
         MRSH  = NBN(IRGN)
         RSHMN = AM12 + RSHRGN(IRGN-1)
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
