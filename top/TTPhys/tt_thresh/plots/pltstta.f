C W/ beam spread and beamstrahlung
CINCLUDE SGTTEF
CINCLUDE UDSRCH
C>>>
CINCLUDE RSHDIS
C*INCLUDE RSHDISN
C>>>
C W/o beam spread and beamstrahlung
C*INCLUDE SGTTEF0
C--
C For perturbative EW.
C*INCLUDE SGTTHR2
C*INCLUDE FHIGGSO
C--
C For pot + perturbative EW.
CINCLUDE SGTTHR2
CINCLUDE FHIGGS
CINCLUDE POTHIG
C--
C For potential approx.
C*INCLUDE SGTTHR2O
C*INCLUDE POTHIG
C
CINCLUDE RKUTTA2
C--
CINCLUDE POTQCD
C*INCLUDE POTRCH
C*INCLUDE POTCPL
C*INCLUDE POTMTN
C>>>
CINCLUDE EI
CINCLUDE F1
CINCLUDE FX
CINCLUDE FH
CINCLUDE GTGAMT
C>>>
C*INCLUDE SETPRM
C For o(alpha_s) correction
CINCLUDE SETPRMC
CINCLUDE TBWCORR
C>>>
C*
C* Main program to draw the ttbar threshold shape.
C* (Update Record)
C*    8/22/91  K.Fujii       Original version for alpha_s dependence.
C*
C*
      IMPLICIT REAL*8 ( A-H, O-Z )
C>>>
C     PARAMETER     ( NPT = 50 )
      PARAMETER     ( NPT = 30 )
C     PARAMETER     ( NPT = 200 )
C>>>
      REAL*8        RDATA(2,0:NPT,0:10)
      CHARACTER*8   JOIN(0:3)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOUNIT   / 21 /
C
C========< Entry Point >================================================
C
C--
C  Set parameters.
C--
C     AMT    = 150.D0
C     AMT    = 170.D0
      AMT    = 175.D0
      VTB2   = 1
      AMH    = 1.D10
      BTH    = 1
C--
C  Alpha_s range.
C--
      NA     = 2
      ALFSMN = 0.11D0
      ALFSMX = 0.13D0
      DA     = (ALFSMX-ALFSMN)/NA
C--
      NA     = 0
      ALFSMN = 0.12D0
C--
C  Energy range.
C--
C>>>
C     RSMN = 2*AMT - 4
C     RSMX = 2*AMT + 1
      RSMN = 2*AMT - 8
      RSMX = 2*AMT + 4
C     RSMN = 2*AMT - 10
C     RSMX = 2*AMT + 4
C>>>
      DRS  = (RSMX-RSMN)/NPT
C--
C     CALL ISRINT(0,2*AMT)
C--
C  Production cross section.
C--
      DO 2000 IA = 0, NA
         MODE = 0
         ALFS = ALFSMN + IA*DA
         DO 200 I = 0, NPT
            RS = RSMN + DRS*I
            CALL SGTTEF(MODE,AMT,ALFS,VTB2,AMH,BTH,RS,SGEFF)
            RDATA(1,I,IA) = RS
            RDATA(2,I,IA) = SGEFF
C>>>
            PRINT *, ' RS, SGEFF = ', RS, SGEFF
C>>>
            MODE       = 3
200      CONTINUE
2000  CONTINUE
C--
C  Prepare TOPDRAW files.
C--
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''( RQQ'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''NEWFRAME'')')
      WRITE(LOUNIT,'(''SET INTENSITY 4'')')
      WRITE(LOUNIT,'(''SET FONT DUPLEX'')')
      WRITE(LOUNIT,'(''SET WINDOW X  2. 12. Y 1.8  9.2'')')
      WRITE(LOUNIT,'(''SET LIMITS X '',2F9.0,'' Y 0 1.5'')') RSMN, RSMX
      WRITE(LOUNIT,'(''SET LABELS SIZE 3.0'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''TITLE 0.5 5.0 ANGLE 90 SIZE 4 ''''S0tt1(pb)'')')
      WRITE(LOUNIT,'(''CASE                          ''''GX  X    '')')
      WRITE(LOUNIT,'(''TITLE 6.0 0.8 SIZE 4 ''''2s0O (GeV)'')')
      WRITE(LOUNIT,'(''CASE                       ''''M UD'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''SET ORDER X Y'')')
C--
      DO 3000 IA = 0, NA
         ALFS = ALFSMN + IA*DA
         WRITE(LOUNIT,'(''( alpha_s = '',1F8.3)') ALFS
         DO 300 I = 0, NPT
            WRITE(LOUNIT,'(5X,5F12.5)') RDATA(1,I,IA), RDATA(2,I,IA)
300      CONTINUE
         WRITE(LOUNIT,'(''JOIN '',A8)') JOIN(IA)
3000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
