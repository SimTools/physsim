C*
C* Main program to draw the ttbar threshold shape.
C* (Update Record)
C*    8/21/91  K.Fujii       Original version.
C*                           W/o radiative correction.
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      PARAMETER     ( NPT = 50 )
      REAL*8        RDATA(2,0:NPT,0:10)
      CHARACTER*8   JOIN(0:4)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT', 'SOLID' /
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
      AMT  = 150.D0
      VFF  = 1.D0
      AMH  = 100.D0
      ALFS = 0.12D0
C--
C  Beta_H range.
C--
      NBT2   = 4
      BTH2MN = 0.5
      BTH2MX = 2.5
      DBT2   = (BTH2MX-BTH2MN)/NBT2
C--
C  Energy range.
C--
      RSMN = 2*AMT - 8
      RSMX = 2*AMT + 4
C     RSMN = 2*AMT - 4
C     RSMX = 2*AMT + 1
      DRS  = (RSMX-RSMN)/NPT
C--
C  Production cross section.
C--
      MODE = 0
      DO 2000 IBT2 = 0, NBT2
         BTH2 = BTH2MN + IBT2*DBT2
         BTH  = SQRT(BTH2)
         IF ( IBT2.EQ.NBT2 ) AMH = 1.D10
         DO 200 I = 0, NPT
            RS = RSMN + DRS*I
            CALL SGTTHR(MODE,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,AMT,VFF,
     .                  AMH,BTH,   RS,SG)
            SG = SG
            RDATA(1,I,IBT2) = RS
            RDATA(2,I,IBT2) = SG
C>>>
C        PRINT *, ' RS, SG = ', RS, SG
C>>>
            MODE = 2
200      CONTINUE
         MODE = 1
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
      WRITE(LOUNIT,'(''SET LIMITS X '',2F9.0,'' Y 0 2.5'')') RSMN, RSMX
      WRITE(LOUNIT,'(''SET TITLE SIZE 3'')')
      WRITE(LOUNIT,'(''TITLE 0.5 5.0 ANGLE 90 SIZE 4 ''''S0tt1(pb)'')')
      WRITE(LOUNIT,'(''CASE                          ''''GX  X    '')')
      WRITE(LOUNIT,'(''TITLE 6.0 0.8 SIZE 4 ''''2s0O (GeV)'')')
      WRITE(LOUNIT,'(''CASE                       ''''M UD'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOUNIT,'(''SET ORDER X Y'')')
C--
      DO 3000 IBT2 = 0, NBT2
         BTH2 = BTH2MN + IBT2*DBT2
         WRITE(LOUNIT,'(''( beta_H^2 = '',F12.5)') BTH2
         DO 300 I = 0, NPT
            WRITE(LOUNIT,'(5X,5F12.5)') RDATA(1,I,IBT2), RDATA(2,I,IBT2)
300      CONTINUE
         WRITE(LOUNIT,'(''JOIN '',A8)') JOIN(IBT2)
3000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
