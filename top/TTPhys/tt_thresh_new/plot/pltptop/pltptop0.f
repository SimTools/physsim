C  C*
C  C* (Update Record)
C  C*   92/03/19  K.Fujii       Parametrization of Cho's new potential
C  C*                           parameters.
C  C*   92/03/19  K.Fujii       This version uses the exact values
C  C*                           obtained by Cho.
C  C*
C
C        SUBROUTINE GTPPAR(X,PPAR)
C
C        IMPLICIT REAL*8 ( A-H, O-Z )
C        REAL*8   X, PPAR(3)
C  C--
C        PARAMETER      ( NPT = 5 )
C        REAL*4           YDAT(NPT,3), DYDAT(NPT,3)
C        DATA YDAT  / 0.20823, 0.23347, 0.23530, 0.20681, 0.16754,
C       .             3.9620 , 3.8077 , 3.7397 , 3.7453 , 3.7485 ,
C       .             0.35914, 0.35398, 0.35719, 0.37349, 0.40100 /
C        DATA DYDAT / 0.01028, 0.01095, 0.00954, 0.00572, 0.00203,
C       .             0.27724, 0.28706, 0.35451, 0.40334, 0.29333,
C       .             0.02136, 0.02130, 0.02466, 0.02769, 0.02078 /
C  C
C  C==================< Entry Point >===================================
C  C
C  C--
C  C  Calculate PPAR.
C  C--
C        IF ( X.GE.0.099D0 .AND. X.LE.0.101D0 ) THEN
C           J = 1
C        ELSE IF ( X.GE.0.109D0 .AND. X.LE.0.111D0 ) THEN
C           J = 2
C        ELSE IF ( X.GE.0.119D0 .AND. X.LE.0.121D0 ) THEN
C           J = 3
C        ELSE IF ( X.GE.0.129D0 .AND. X.LE.0.131D0 ) THEN
C           J = 4
C        ELSE IF ( X.GE.0.139D0 .AND. X.LE.0.141D0 ) THEN
C           J = 5
C        ELSE
C           PRINT *, ' ERROR in GTPPAR: ALFS = ', X
C        ENDIF
C  C>>>
C  C     PPAR(1) = YDAT(J,1) - 5*DYDAT(J,1)
C        PPAR(1) = YDAT(J,1) + 5*DYDAT(J,1)
C  C     PPAR(1) = YDAT(J,1)
C  C>>>
C        PPAR(2) = YDAT(J,2)
C        PPAR(3) = YDAT(J,3)
C  C--
C  C  That's it.
C  C--
C        RETURN
C        END
C
C=======================================================================
C  Subroutines necessary to run this job.
C=======================================================================
C
C--
C  QCD corrections.
C--
CINCLUDE GTCEFF
CINCLUDE GMTWMB
C>>>
C*INCLUDE QCDPOT1
C>>>
CINCLUDE QCDPOT
CINCLUDE ADJ1
CINCLUDE ADJ2
CINCLUDE QCDS1
CINCLUDE QCDS2
CINCLUDE QCDS3
CINCLUDE BNDCND
CINCLUDE RKUTAS
CINCLUDE QCDP1
CINCLUDE QCDP2
CINCLUDE QCDP3
CINCLUDE RKUTAP
CINCLUDE DELFN
CINCLUDE RUNWID
CINCLUDE LNT0FI
CINCLUDE LNT1FI
CINCLUDE NTGRD1
CINCLUDE NTGRD2
CINCLUDE NTGRD3
CINCLUDE LNT0WA
CINCLUDE PHS0
CINCLUDE LNT2WA
CINCLUDE PHS2
CINCLUDE TBWCORR
CINCLUDE FX
CINCLUDE GTQCDP
CINCLUDE GAMMA0
C--
C  Cross section and asymmetry.
C--
CINCLUDE DSGDP0
CINCLUDE GRQDSG
CINCLUDE GRQINT
CINCLUDE GRQADJ
CINCLUDE GRQPTP
C
C=======================================================================
C  Main programs.
C=======================================================================
C
C*
C*  93/04/19  K.Fujii        This program plots ds/dp.
C*
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
C--
      COMMON /KFLAGS/ KFNOFI 
      INTEGER*4       KFNOFI
C--
      PARAMETER     ( MXxHPK = 40000 )
      COMMON /PAWC/ M(MXxHPK)
      DATA LOU    /21/
C
C========< Entry Point >================================================
C
C--
C  Set constants.
C--
      ALP  = 1/128.D0
      SN2W = 0.23D0
      AMSZ = 91.17D0
      AMSW = 80.0D0
      AMSB =  5.0D0
      AMSH =  1.D10
      BETH =  0
C--
C  Switch off FI corrections if you want.
C     KFNOFI = (0,1) = (ON,OFF)
C--
C      KFNOFI = 1
C--
C  Set top parameters.
C--
      ALPS = 0.12
C      AMST = 150.D0
      AMST = 170.D0
      VTB2 = 1.D0
C--
C  Energy range.
C--
C>>>
      NE   = 40
      EMN  = -6.0
      EMX  = +2.0
      DE   = (EMX-EMN)/NE
C>>>
C     EMN  = -2.5514D0
C     EMN  = 0.D0
C     EMN  = +1.D0
C
C DELE = +2 for mt=150(new running width)
C      EMN  = -0.5514D0
C 94/08/03
C DELE = +2 for mt=170(new running width)
      EMN  = -0.595
C>>>
      NE   = 0
C--
C  Momentum range.
C--
      NP   = 50
      PMN  = 0.D0
      PMX  = 50.D0
      DP   = (PMX-PMN)/NP
C--
C  Prepare TOPDRAW files.
C--
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''( dsig/dp'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''NEWFRAME'')')
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''SET INTENSITY 4'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET WINDOW X  2.5 12.5 Y 1.8  9.2'')')
      WRITE(LOU,'(''SET LIMIT X '',2E15.8)') PMN, PMX
      WRITE(LOU,'(''SET LIMIT Y 0. 0.125'')')
      WRITE(LOU,'(''(SET SCALE Y LOG'')')
      WRITE(LOU,'(''SET LABELS SIZE 3.0'')')
C--
      WRITE(LOU,'(''TITLE 0.5 3.2 ANGLE 90 SIZE 4 ''''dS0tt1/d|p|'',
     .            '' (pb/GeV)'')')
      WRITE(LOU,'(''CASE                          '''' GX  X     '')')
      WRITE(LOU,'(''TITLE 6.5 0.8 SIZE 4 ''''|p| (GeV)'')')
      WRITE(LOU,'(''('')')
C--
      WRITE(LOU,'(''( AMST = '',1F7.3)') AMST
      WRITE(LOU,'(''( ALFS = '',1F7.3)') ALPS
      WRITE(LOU,'(''( AMSH = '',1E12.3)') AMSH
      WRITE(LOU,'(''( BETH = '',1F7.3)') BETH
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET ORDER DUMMY X Y'')')
C--
C  Plot ds/dp.
C--
      MODE = 0
      DO 1000 IE = 0, NE
         E = EMN + DE*IE
         SGT  = 0
         DO 100 IP = 0, NP
            IF ( IP.EQ.0 .OR. IP.EQ.NP ) THEN
               WT = 1
            ELSE IF ( MOD(IP,2).EQ.0 ) THEN
               WT = 2
            ELSE
               WT = 4
            ENDIF
            P = PMN + IP*DP
            CALL DSGDP0(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,
     .                  E,P,DSG)
            SGT = SGT + DSG*WT
            WRITE(LOU,*) E, P, DSG
            PRINT *, 'P = ', REAL(P), ' ds/dp = ', REAL(DSG)
            MODE = 4
100      CONTINUE
         SGT = SGT*DP/3
         PRINT *, 'E = ', E, ' SGT = ', SGT
         WRITE(LOU,'(''JOIN '')')
         MODE = 3
1000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
