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
C  C     PPAR(1) = YDAT(J,1)
C        PPAR(1) = YDAT(J,1) + 5*DYDAT(J,1)
C  C>>>
C        PPAR(2) = YDAT(J,2)
C        PPAR(3) = YDAT(J,3)
C  C--
C  C  That's it.
C  C--
C        RETURN
C        END
CINCLUDE SGTTHR2
CINCLUDE RKUTTA2
CINCLUDE EI
CINCLUDE FX
CINCLUDE FH
CINCLUDE F1
CINCLUDE GTGAMT
CINCLUDE FHIGGS
C*INCLUDE FHIGGSO
C>>>
C*INCLUDE SETPRM
C for o(alfs) correction
CINCLUDE SETPRMC
CINCLUDE TBWCORR
C>>>
CINCLUDE POTHIG
C*INCLUDE POTRCH
C*INCLUDE QCDPOT1
CINCLUDE POTQCD
C*INCLUDE POTQCDO
C*INCLUDE POTCPL
C*INCLUDE POTMTN
C>>>
C*
C* Main program to draw the ttbar threshold shape.
C* (Update Record)
C*    8/21/91  K.Fujii       Original version.
C*                           W/o radiative correction.
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
C     PARAMETER     ( NPT = 200 )
      PARAMETER     ( NPT = 60 )
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
      ALF  = 1.D0/128
      S2W  = 0.230D0
      AMZ  = 91.17D0
      GMZ  = 2.5D0
      AMW  = 80.D0
      AMB  = 5.D0
C--
C     AMT  = 150.D0
      AMT  = 170.D0
      VFF  = 1.D0
C     VFF  = SQRT(0.1D0)
      AMH  = 1.D10
      BTH  = 1.D0
C--
C  Alpha_s range.
C--
      NA     = 2
      ALFSMN = 0.11D0
      ALFSMX = 0.13D0
      DA     = (ALFSMX-ALFSMN)/NA
C>>>
      NA     = 0
      ALFSMN = 0.12D0
C>>>
C--
C  Energy range.
C--
C     RSMN = 2*AMT - 10
C     RSMX = 2*AMT + 4
      RSMN = 2*AMT - 8
      RSMX = 2*AMT + 4
C     RSMN = 2*AMT - 4
C     RSMX = 2*AMT + 1
      DRS  = (RSMX-RSMN)/NPT
C--
      EMN  = RSMN - 2*AMT
      EMX  = RSMX - 2*AMT
C--
C  Production cross section.
C--
      DO 2000 IA = 0, NA
         ALFS = ALFSMN + IA*DA
         MODE = 0
         DO 200 I = 0, NPT
            RS = RSMN + DRS*I
            CALL SGTTHR(MODE,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,AMT,VFF,
     .                  AMH,BTH,   RS,SG)
            SG = SG
C>>>
C           RDATA(1,I,IA) = RS
            RDATA(1,I,IA) = RS - 2*AMT
C>>>
            RDATA(2,I,IA) = SG
C>>>
C        PRINT *, ' RS, SG = ', RS, SG
C>>>
            MODE       = 2
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
C     WRITE(LOUNIT,'(''SET LIMITS X '',2F9.0,'' Y 0 2.5'')') RSMN, RSMX
      WRITE(LOUNIT,'(''SET LIMITS X '',2F9.0,'' Y 0 2.5'')') EMN, EMX
      WRITE(LOUNIT,'(''SET TITLE SIZE 3'')')
      WRITE(LOUNIT,'(''TITLE 0.5 5.0 ANGLE 90 SIZE 4 ''''S0tt1(pb)'')')
      WRITE(LOUNIT,'(''CASE                          ''''GX  X    '')')
C     WRITE(LOUNIT,'(''TITLE 6.0 0.8 SIZE 4 ''''2s0O (GeV)'')')
C     WRITE(LOUNIT,'(''CASE                       ''''M UD'')')
      WRITE(LOUNIT,'(''TITLE 6.0 0.8 SIZE 4 ''''E (GeV)'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOUNIT,'(''SET ORDER X Y'')')
C--
      DO 3000 IA = 0, NA
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
