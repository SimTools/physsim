C*
C* Main program to draw the ttbar threshold shape.
C* (Update Record)
C*    8/22/91  K.Fujii       Original version for M_H dependence.
C*
C*
      IMPLICIT REAL*8 ( A-H, O-Z )
      PARAMETER     ( NPT = 50 )
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
      AMT    = 170.D0
      ALFS   = 0.12D0
      VTB2   = 1
      BTH    = 1
C--
C  M_H range.
C--
      NA     = 2
      ALFSMN = 50.D0
      ALFSMX = 150.D0
      DA     = (ALFSMX-ALFSMN)/NA
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
C  Production cross section.
C--
      DO 2000 IA = 0, NA
         MODE = 0
         AMH  = ALFSMN + IA*DA
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
         AMH  = ALFSMN + IA*DA
         WRITE(LOUNIT,'(''( M_H = '',1F8.3)') AMH
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
