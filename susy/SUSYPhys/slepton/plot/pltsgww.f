CX*INCLUDE RSHDIS
CINCLUDE RSHDISN
CINCLUDE UDSRCH
C--
C*
C*
C*  Total cross section for e+e- --> X+ X-.
C*    RS = sqrt(S)
C*
C*  5/16/91  K.Fujii       This version can handle X+1 X-2 case.
C*
      IMPLICIT    REAL*8  ( A-H, O-Z )
C--
      CHARACTER*7 JOIN(1:4)
      DATA JOIN   / 'SOLID  ',
     .              'DOTDASH',
     .              'DASH   ',
     .              'DOT    '/
      DATA LOU / 20 /
C
C========< Entry Point >================================================
C
C--
C  Total cross section.
C--
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET WINDOW X 2.5 12.5 Y 2.2 9.6'')')
      WRITE(LOU,'(''SET LIMITS X 100. 600. Y 0.01 100.'')')
      WRITE(LOU,'(''SET SCALE Y LOGARITHMIC'')')
      WRITE(LOU,'(''SET LABELS SIZE 3.0'')')
      WRITE(LOU,'(''SET TITLE  SIZE 3.0'')')
      WRITE(LOU,'(''TITLE   7. 8.8 SIZE 3.5'',
     .            '' '''' e2+3e2-3 --> W2+3W2-3'')')
      WRITE(LOU,'(''CASE                   '',
     .            '' ''''  X X X X      X X X X'')')
      WRITE(LOU,'(''TITLE 6.8 1.0 SIZE 4 ''''2s0O (GeV)'')')
      WRITE(LOU,'(''CASE          ''''M UD      '')')
      WRITE(LOU,'(''TITLE 0.5 5.0 ANGLE 90. SIZE 4 '',
     .                            '' ''''S (pb)'')')
      WRITE(LOU,'(''CASE             ''''G     '')')
      WRITE(LOU,'(''SET ORDER X DUMMY Y'')')
C--
C  sqrt(s) range.
C--
C     NRS  = 100
C     RSMN = 160.1
C--
      NRS  = 0
      RSMN = 350.
C--
      DRS  =  5
      MODE = 0
      DO 10 IRS = 0, NRS
         RS   = RSMN + IRS*DRS
         CALL SGWWEF(MODE,RS,SG)
         CALL SGTWWA(RS,SG0)
         WRITE(LOU,'(3E15.5)') RS, SG0, SG
         PRINT *, ' RS = ', RS, ' GeV  sig0 = ', SG0, ' pb',
     .                          ' sig = ', SG, ' pb'
         MODE = 3
10    CONTINUE
      WRITE(LOU,'(''JOIN '',A7)') JOIN(1)
C--
      STOP
      END
C*
C* (Update Record)
C*   92/06/26  K.Fujii       Reduced # bins for sig_0 tabulations.
C*
 
      SUBROUTINE SGWWEF(MODE,RS,SGEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE
      REAL   *8    RS, SGEFF
C--
      REAL   *4    XRS
C--
      PARAMETER    ( NRSH = 200 )
      INTEGER*4    NBN(3)
      DATA NBN     /  50, 20, 100 /
      REAL   *8    RSHRGN(0:3), DRSHRG(3), SGDAT(0:NRSH,3)
      DATA RSHRGN  /   0.D0,  5.D0, 25.D0, 1025.D0 /
      DATA NCALL   /0/
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         NBM    = 100
         DBM    = 1./NBM
         AMW    = 80
         AMX12  = 2*AMW
      ENDIF
C--
C  Check if RS is in the acceptable range.
C--
      E = RS - 2*AMW
      IF ( E.LT.RSHRGN(0) .OR. E.GT.RSHRGN(3) ) THEN
         PRINT *, ' SGWWEF: sqrt(s) = ', RS, ' is out of range.'
         PRINT *, ' AMX = ', AMW, ' EMN = ', RSHRGN(0),
     .            ' EMX = ', RSHRGN(3)
         PRINT *, '       : causes forced STOP.'
         STOP
      ENDIF
C--
C  Tabulate cross section w/o RC.
C--
      IF ( MODE.LE.2 ) THEN
         PRINT *, 'SGWWEF now tabulates sigma_0.'
C--
         DO 100 IRGN = 1, 3
            MRSH  = NBN(IRGN)
            RSHMN = AMX12 + RSHRGN(IRGN-1)
            RSHMX = AMX12 + RSHRGN(IRGN)
            DRSH  = (RSHMX-RSHMN)/MRSH
            DRSHRG(IRGN) = DRSH
            DO 10 IRSH = 0, MRSH
               RSH = RSHMN + DRSH*IRSH
               CALL SGTWWA(RSH,SG)
               SGDAT(IRSH,IRGN) = SG
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
C*
C*  Integrated cross section for e+e- --> WW.
C*    RS = sqrt(S)
C*
      SUBROUTINE SGTWWA(RS,SGT)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
      REAL*8   RS, SGT
C
C========< Entry Point >================================================
C
C--
C  Integrate differential cross section.
C--
      NY   = 1000
      YMN  = -1
      YMX  =  1
      DY   = (YMX-YMN)/NY
      SGT = 0
      DO 10 IY = 0, NY
         Y = YMN + IY*DY
         IF ( IY.EQ.0 .OR. IY.EQ.NY ) THEN
            WT = 1
         ELSE IF ( MOD(IY,2).EQ.0 ) THEN
            WT = 2
         ELSE
            WT = 4
         ENDIF
         CALL SGWWNP(RS,Y,SG)
         SGT = SGT + WT*SG
10    CONTINUE
      SGT = SGT*DY/3
C--
C  That's it.
C--
      RETURN
      END
 
      SUBROUTINE SGWWNP(RS,Y,SG)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      PARAMETER  ( IxPRC = 2 )
      REAL   *8  RS, Y, SG
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xPI    = ACOS(-1.D0)
         x2PI   = 2*xPI
         x4PI   = 4*xPI
         xGV2PB = 3.893796623D8
C        xALF   = 1.D0/137.035989561D0
         xALF   = 1.D0/128.D0
C--
         AMZ    = 91.17D0
         GMZ    =  2.5D0
C        GMZ    =  .0D0
         AMW    = 80.0D0
         GMW    =  2.0D0
         AMZ2   = AMZ*AMZ
         AMW2   = AMW*AMW
         AMZ4   = AMZ2*AMZ2
C--
C        xSIN2W = 0.228842366686D0
         xSIN2W = 1 - AMW2/AMZ2
         xCOS2W = 1 - xSIN2W
         xSINW  = SQRT(xSIN2W)
         xCOSW  = SQRT(xCOS2W)
         xCOTW  = xCOSW/xSINW
C--
         C0     = SQRT(x4PI*xALF)
         CV     = C0/(xSINW*xCOSW)*(-0.25D0+xSIN2W)
         CA     = C0/(xSINW*xCOSW)*( 0.25D0 )
         CWWZ   = C0*xCOTW
         CW     = C0/xSINW/2/SQRT(2.D0)
         FAC    =  xGV2PB/(64*xPI**2)
      ENDIF
C--
C  Set independent variables.
C--
      EBM  = RS/2
      BT   = (EBM-AMW)*(EBM+AMW)
      IF ( BT.LE.0.D0 ) THEN
         SG = 0
         RETURN
      ELSE
         BT   = SQRT(BT)/EBM
      ENDIF
      S    = RS*RS
      T    = AMW2 - (S/2)*(1-BT*Y)
      P1P2 = S/2
      P1Q1 = (AMW2-T)/2
C--
C  Set independent variables.
C--
      CALL       TTWWNP(AMZ,GMZ,AMW,GMW,
     .                  C0,CV,CA,CWWZ,CW,
     .                  P1P2,P1Q1,TTA)
      SG = FAC*TTA*BT/S
C--   ds/dcos(theta)
      SG = SG*2*xPI
C--
C  That's it.
C--
      RETURN
      END
 
      SUBROUTINE TTWWNP(AMZ,GMZ,AMW,GMW,
     .                  C0,CV,CA,CWWZ,CW,
     .                  P1P2,P1Q1,TTA)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
C--
C  Statement functions.
C--
      RBW(A,B) =  A-B*B
      BW2(A,B) = (A-B*B)**2+B*B*GMZ*GMZ
C--
C  Ddifferential cross section is given by the follwing.
C     d(sigma)/d(omega) = TTA*(e**2/(sin2w*cos2w))**2*BT/(64*xPI**2*S)
C                       = TTA*(xALF**2/(sin2w*cos2w))**2*BT/(4*S)
C  REDUCE output here.
C--
      TTAGG=C0**4*(-3.*AMW**6*P1P2-6.*AMW**4*P1Q1**2+6.*AMW**4*
     .      P1Q1*P1P2-6.*AMW**4*P1P2**2+4.*AMW**2*P1Q1**2*P1P2-4.*AMW
     .      **2*P1Q1*P1P2**2+3.*AMW**2*P1P2**3-2.*P1Q1**2*P1P2**2+2.*
     .      P1Q1*P1P2**3)/(AMW**4*P1P2**2)
      TTATT=8.*CW**4*(-AMW**6*P1P2-10.*AMW**4*P1Q1**2-2.*AMW**4*
     .      P1Q1*P1P2+8.*AMW**2*P1Q1**3+4.*AMW**2*P1Q1**2*P1P2-8.*
     .      P1Q1**4+8.*P1Q1**3*P1P2)/(AMW**8-4.*AMW**6*P1Q1+4.*AMW**4
     .      *P1Q1**2)
      TTAZG=4.*C0**2*CV*CWWZ*(3.*AMW**6*P1P2+6.*AMW**4*P1Q1**2-
     .      6.*AMW**4*P1Q1*P1P2+6.*AMW**4*P1P2**2-4.*AMW**2*P1Q1**2*
     .      P1P2+4.*AMW**2*P1Q1*P1P2**2-3.*AMW**2*P1P2**3+2.*P1Q1**2*
     .      P1P2**2-2.*P1Q1*P1P2**3)/(2.*AMW**4*P1P2**2-AMW**4*P1P2*
     .      AMZ**2)
      TTAGT=4.*C0**2*CW**2*(-3.*AMW**6*P1P2-6.*AMW**4*P1Q1**2-3.
     .      *AMW**4*P1P2**2+4.*AMW**2*P1Q1**3-2.*AMW**2*P1Q1**2*P1P2+
     .      4.*AMW**2*P1Q1*P1P2**2-4.*P1Q1**3*P1P2+4.*P1Q1**2*P1P2**2)
     .      /(AMW**6*P1P2-2.*AMW**4*P1Q1*P1P2)
      TTAZZ=4.*(CV**2+CA**2)*CWWZ**2*(-3.*AMW**6*P1P2-6.*AMW**4*P1Q1**2+
     .      6.*AMW**4*P1Q1*P1P2-6.*AMW**4*P1P2**2+4.*AMW**2*P1Q1**2*
     .      P1P2-4.*AMW**2*P1Q1*P1P2**2+3.*AMW**2*P1P2**3-2.*P1Q1**2*
     .      P1P2**2+2.*P1Q1*P1P2**3)/(4.*AMW**4*P1P2**2-4.*AMW**4*
     .      P1P2*AMZ**2+AMW**4*AMZ**4)
      TTAZT=8.*(CV-CA)*CW**2*CWWZ*(3.*AMW**6*P1P2+6.*AMW**4*P1Q1**2+
     .      3.*AMW**4*P1P2**2-4.*AMW**2*P1Q1**3+2.*AMW**2*P1Q1**2*P1P2
     .      -4.*AMW**2*P1Q1*P1P2**2+4.*P1Q1**3*P1P2-4.*P1Q1**2*P1P2**2)
     .      /(2.*AMW**6*P1P2-AMW**6*AMZ**2-4.*AMW**4*P1Q1*P1P2+2.*
     .      AMW**4*P1Q1*AMZ**2)
      TTA = TTAZZ + TTAGG + TTATT + TTAZG + TTAZT + TTAGT
C--
C  That's it.
C--
      RETURN
      END
