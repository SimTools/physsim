C*
C* Main program to prepare sig_eff table.
C* (Update Record)
C*    8/22/91  K.Fujii       Original version for !V_tb!^2 dependence.
C*
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      DATA LOU   / 20 /
C
C========< Entry Point >================================================
C
C--
C  Allocate output data set.
C--
C      CALL UALCPS(LOU,'TKSF.@.@SGTBL@.MHBT.DATA','RENEW','FB',30,IRT)
C--
C  Set parameters.
C--
C     AMT    = 150.D0
      AMT    = 170.D0
      ALFS   = 0.12D0
      VTB2   = 1.D0
C--
C  Beta_H^2 range.
C--
      NAL    = 6
      ALMN   = 1.D-3
      ALMX   = 4.D0
      DAL    = (ALMX-ALMN)/NAL
C--
C  M_H range.
C--
      NAM    = 10
      AMMN   = 50.D0
      AMMX   = 300.D0
      DAM    = (AMMX-AMMN)/NAM
C--
C  Energy range.
C--
      NRS    = 10
      RSMN   = 2*AMT - 7
      RSMX   = 2*AMT + 3
      DRS    = (RSMX-RSMN)/NRS
C--
C  Production cross section.
C--
      MODE = 0
      DO 1000 IAL = 0, NAL
         BTH2 = ALMN + IAL*DAL
         BTH  = SQRT(BTH2)
         DO 100 IAM = 0, NAM
            AMH = AMMN + IAM*DAM
            DO 10 I = 0, NRS
               RS = RSMN + DRS*I
               CALL SGTTEF(MODE,AMT,ALFS,VTB2,AMH,BTH,RS,SGEFF)
               SGEFF = SGEFF
               PRINT *, BTH2, AMH, RS, SGEFF
               WRITE(LOU,'(4E15.6)') BTH2, AMH, RS, SGEFF
               MODE       = 3
10          CONTINUE
            MODE = 1
100      CONTINUE
1000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
