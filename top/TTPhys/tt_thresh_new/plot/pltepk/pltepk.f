C*
C* 93/04/19  K.Fujii                Plots E_peak vs m_t.
C*
C*
      IMPLICIT REAL*8 ( A-H, O-Z )
      DATA LOU /21/
 
      ALF    = 1/128.0D0
      S2W    =  0.230D0
      AMB    =  5.0D0
      AMW    = 80.0D0
      AMZ    = 91.17D0
C--
      VTB2 = 1
      AMH  = 1.D10
      BTH  = 0.D10
      ALFS = 0.12D0
      EPK  = -2.8D0
C--
      NM    = 40
      AMTMN = 120
      AMTMX = 180
      DM    = (AMTMX-AMTMN)/NM
C--
      MODE = 0
C--
      WRITE(LOU,*) 'SET ORDER DUMMY X Y DUMMY'
      DO 100 IM = 0, NM
         AMT = AMTMN + IM*DM
         CALL EPKHNT( MODE, ALF,S2W,AMZ,AMW,AMB,
     .                      AMT,ALFS,VTB2,AMH,BTH, EPK,SGPK, NTRY )
         IF ( NTRY.LT.0 ) THEN
            PRINT *, '>>>> EPKHNT failed to find peak >>>>'
            PRINT *, '   AMT   = ', AMT
            PRINT *, '   ALFS  = ', ALFS
            PRINT *, '   VTB2  = ', VTB2
            PRINT *, '   AMH   = ', AMH
            PRINT *, '   BTH   = ', BTH
         ENDIF
         WRITE(LOU,*) REAL(ALFS), REAL(AMT), REAL(EPK), REAL(SGPK)
C>>>
         PRINT *, ' NTRY, AMT, EPK, SGPK = ',
     .              NTRY, REAL(AMT), REAL(EPK), REAL(SGPK)
C>>>
         MODE = 1
100   CONTINUE
      WRITE(LOU,*) 'JOIN SOLID'
      END
