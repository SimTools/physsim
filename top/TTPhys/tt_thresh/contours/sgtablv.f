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
C      CALL UALCPS(LOU,'TKSF.@.@SGTBL@.MTV2.DATA','RENEW','FB',30,IRT)
C--
C  Set parameters.
C--
C     AMT    = 150.D0
C     AMT    = 170.D0
C     AMT    = 175.D0
      AMT    = 165.D0
      ALFS   = 0.12D0
      AMH    = 1.D10
      AMZ    = 91.1876D0
      BTH    = 1
C--
C  MMODE = (0,1) = (onshell,msbar)
C--
      MMODE = 1
C--
C  !V_tb!^2 range.
C--
C>>>
C     NAL    = 10
C     ALMN   = 0.510D0
C     ALMX   = 1.490D0
      NAL    = 14
      ALMN   = 0.300D0
      ALMX   = 1.700D0
C>>>
      DAL    = (ALMX-ALMN)/NAL
C--
C  m_t range.
C--
C>>>
C     NAM    = 20
C     AMMN   = AMT - 0.2D0
C     AMMX   = AMT + 0.2D0
C     NAM    = 40
C     AMMN   = AMT - 0.4D0
C     AMMX   = AMT + 0.4D0
      NAM    = 40
      AMMN   = AMT - 0.2D0
      AMMX   = AMT + 0.2D0
C>>>
      DAM    = (AMMX-AMMN)/NAM
C--
C  Energy range.
C--
      AMTP   = AMTPOL(MMODE,AMT,ALFS,AMZ)
      NRS    = 10
      RSMN   = 2*AMTP - 7
      RSMX   = 2*AMTP + 3
      DRS    = (RSMX-RSMN)/NRS
C--
C  Production cross section.
C--
      MODE = 0
      DO 1000 IAL = 0, NAL
         VTB2 = ALMN + IAL*DAL
         DO 100 IAM = 0, NAM
            AMT  = AMMN + IAM*DAM
            AMTP = AMTPOL(MMODE,AMT,ALFS,AMZ)
            DO 10 I = 0, NRS
               RS = RSMN + DRS*I
               CALL SGTTEF(MODE,AMTP,ALFS,VTB2,AMH,BTH,RS,SGEFF)
               SGEFF = SGEFF
               PRINT *, VTB2, AMTP, RS, SGEFF
               WRITE(LOU,'(4E17.8)') VTB2, AMTP, RS, SGEFF
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
