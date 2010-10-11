CINCLUDE COUPLSUB
CINCLUDE SFMASS
CINCLUDE INOMIX
CINCLUDE SWPELM
CINCLUDE INIPRM
CINCLUDE DSSESE
CINCLUDE SGSESE
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  AM(2), AMX(4), CGE(2), CZE(2), CGSE(2), CZSE(2)
      COMPLEX*16 CX(2,2,4)
      PARAMETER  ( NP = 1000 )
      REAL   *8  XDATA(0:NP,4), YDATA(0:NP,4)
C
C==================< Entry Point >======================================
C
C--
C  Set experimental parameters.
C--
      POLE = +1
      RSMX = 500
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
      TNB   = 2.D0
C--
      x2PI  = 2*ACOS(-1.D0)
C--
      CALL INIPRM(4,IFLG,AM0,AMU,AM2,TNB,
     .            AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE)
C--
      WRITE(20,'(''( AM0  = '',F8.2)') AM0
      WRITE(20,'(''( AMU  = '',F8.2)') AMU
      WRITE(20,'(''( AM2  = '',F8.2)') AM2
      WRITE(20,'(''( TNB  = '',F8.2)') TNB
      WRITE(20,'(''( AMEL = '',F8.2)') AM(1)
      WRITE(20,'(''( AMER = '',F8.2)') AM(2)
      WRITE(20,'(''( AMX1 = '',F8.2)') AMX(1)
      WRITE(20,'(''( AMX2 = '',F8.2)') AMX(2)
      WRITE(20,'(''( AMX3 = '',F8.2)') AMX(3)
      WRITE(20,'(''( AMX4 = '',F8.2)') AMX(4)
C--
C  Loop over final state chirality combinations.
C--
      AMEL = AM(1)
      AMER = AM(2)
      DO 1000 MODE = 1, 4
C--
C  Set selectron masses.
C--
         IF ( MODE.EQ.1 ) THEN
            AM(1) = AMEL
            AM(2) = AMEL
         ELSE IF ( MODE.EQ.2 ) THEN
            AM(1) = AMER
            AM(2) = AMER
         ELSE IF ( MODE.EQ.3 ) THEN
            AM(1) = AMER
            AM(2) = AMEL
         ELSE IF ( MODE.EQ.4 ) THEN
            AM(1) = AMEL
            AM(2) = AMER
         ENDIF
C--
C  Loop over sqrt(s).
C--
         NRS  = 100
         RSMN = AM(1) + AM(2) + 1.D-3
         DRS  = (RSMX-RSMN)/NRS
         CSMX = 1.D0
C--
         DO 100 IRS = 0, NRS
            RS = RSMN + IRS*DRS
            CALL SGSESE(MODE,AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE,
     .                  RS,POLE, CSMX, SGT)
            XDATA(IRS,MODE) = RS
            YDATA(IRS,MODE) = SGT
C           PRINT *, ' RS = ', RS, '(GeV) SGT = ', SGT, '(pb)'
100      CONTINUE
1000  CONTINUE
C--
C  Prepare topdraw data.
C--
      WRITE(20,'(''NEW FRAME'')')
      WRITE(20,'(''SET FONT DUPLEX'')')
      WRITE(20,'(''('')')
      WRITE(20,'(''SET WINDOW X 2.5 12.5 Y 2.2 9.6'')')
      WRITE(20,'(''SET LIMITS X 270 500 Y 5.E-4 5.E+0'')')
      WRITE(20,'(''SET SCALE Y LOG                   '')')
      WRITE(20,'(''SET LABELS SIZE 3      '')')
      WRITE(20,'(''('')')
      WRITE(20,'(''TITLE 6.4 1.2 SIZE 5 ''''2s0O (GeV)'')')
      WRITE(20,'(''CASE                 ''''M UD      '')')
      WRITE(20,'(''TITLE 0.4 4.7 SIZE 5 ANGLE 90 ''''S0tot1 (pb) '')')
      WRITE(20,'(''CASE                          ''''GX   X      '')')
      WRITE(20,'(''('')')
      WRITE(20,'(''SET ORDER X Y '')')
      DO 200 IRS = 0, NRS
         WRITE(20,'(2E15.5)') XDATA(IRS,1), YDATA(IRS,1)
200   CONTINUE
      WRITE(20,'(''JOIN DASH'')')
      DO 300 IRS = 0, NRS
         WRITE(20,'(2E15.5)') XDATA(IRS,2), YDATA(IRS,2)
300   CONTINUE
      WRITE(20,'(''JOIN SOLID'')')
      DO 400 IRS = 0, NRS
         WRITE(20,'(2E15.5)') XDATA(IRS,3), YDATA(IRS,3)
400   CONTINUE
      WRITE(20,'(''JOIN DOTDASH'')')
      DO 500 IRS = 0, NRS
         WRITE(20,'(2E15.5)') XDATA(IRS,3), YDATA(IRS,4)
500   CONTINUE
      WRITE(20,'(''JOIN DOT'')')
C--
C  That's it.
C--
      STOP
      END
