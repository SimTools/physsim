C*
C* Main program to prepare sig_eff table.
C* (Update Record)
C*    8/22/91  K.Fujii       Original version for alpha_s dependence.
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
C      CALL UALCPS(LOU,'TKSF.@.@SGTBL@.MTAL.DATA','RENEW','FB',30,IRT)
C--
C  Set parameters.
C--
C     AMT    = 150.D0
      AMT    = 170.D0
      VTB2   = 1.0
      AMH    = 1.D10
      BTH    = 1
C--
C  Alpha_s range.
C--
      NAL    = 12
      ALMN   = 0.100D0
      ALMX   = 0.136D0
      DAL    = (ALMX-ALMN)/NAL
C--
C  m_t range.
C--
      NAM    = 80
      AMMN   = AMT - 0.5
      AMMX   = AMT + 0.5
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
      DO 1000 IAL = 0, NAL
         ALFS = ALMN + IAL*DAL
         MODE = 0
         DO 100 IAM = 0, NAM
            AMT = AMMN + IAM*DAM
            DO 10 I = 0, NRS
               RS = RSMN + DRS*I
               CALL SGTTEF(MODE,AMT,ALFS,VTB2,AMH,BTH,RS,SGEFF)
               SGEFF = SGEFF
               PRINT *, ALFS, AMT, RS, SGEFF
               WRITE(LOU,'(4E15.6)') ALFS, AMT, RS, SGEFF
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
