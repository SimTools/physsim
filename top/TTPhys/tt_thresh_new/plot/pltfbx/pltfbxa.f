CC**********************************************************************
C*
C* (Purpose)
C*     This program calculates the FB-asymmetry as a function of
C*     alpha_s at a given delta_E value.
C* (Update Record)
C*    93/04/19  K.Fujii                Original version w/ ISR.
C*
CC**********************************************************************

      PROGRAM pltfbxa
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      DATA LOU /21/
C
C========< Entry Point >================================================
C
C--
C  Allocate output file.
C--
C     CALL UALCPS(LOU,'TKSF.@@.@FBDX@.ALFS.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@FBDX@.ALFS.BM.TDR','RENEW','FB',30,IRT)
C--
C  Initialize constants.
C--
      VTB2   = 1.0D0
      AMH    = 1.D10
      BTH    = 0.D10
C--
C>>>
      ALFSMN = 0.10D0
      ALFSMX = 0.14D0
      NA     = 1
      DA     = (ALFSMX-ALFSMN)/NA
C>>>
C--
C     AMT    = 150.D0
      AMT    = 170.D0
C--
      EPK    = -2.5D0
C--
      ALF    = 1/128.D0
      S2W    = 0.23D0
      AMZ    = 91.17D0
      AMW    = 80.0D0
      AMB    =  5.0D0
C--
C  Energy.
C--
      DE     = 1.0D0
C--
C  Initialize ISR.
C--
      CALL ISRINT(0,2*AMT)
C--
C  Topdraw commands.
C--
      WRITE(LOU,*) '( '
      WRITE(LOU,*) '( delta_E  = ', REAL(DE)
      WRITE(LOU,*) '( m_t      = ', REAL(AMT)
      WRITE(LOU,*) '( !V_tb!^2 = ', REAL(VTB2)
      WRITE(LOU,*) '( No Higgs '
      WRITE(LOU,*) '( '
      WRITE(LOU,*) 'SET LIMITS X ', ALFSMN, ALFSMX, ' Y -0.05 0.15'
      WRITE(LOU,*) 'SET ORDER X DUMMY Y'
C--
C  Loop over alpha_s.
C--
      MODE2 = 0
      DO 1000 IA = 0, NA
         ALFS = ALFSMN + IA*DA
         WRITE(LOU,*) '( alpha_s  = ', REAL(ALFS)
         CALL EPKHNT( 0, ALF,S2W,AMZ,AMW,AMB,
     .                   AMT,ALFS,VTB2,AMH,BTH, EPK,SGPK, NTRY )
         IF ( NTRY.LT.0 ) THEN
            PRINT *, '>>>> EPKHNT failed to find peak >>>>'
            PRINT *, '   AMT   = ', AMT
            PRINT *, '   ALFS  = ', ALFS
            PRINT *, '   VTB2  = ', VTB2
            PRINT *, '   AMH   = ', AMH
            PRINT *, '   BTH   = ', BTH
         ELSE
            WRITE(LOU,*) '( E_peak   = ', REAL(EPK)
            WRITE(LOU,*) '( SG0_peak = ', REAL(SGPK)
            E     = EPK + DE
            CALL SGAFBF(MODE2,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                  E,SGT,AFB)
            MODE2 = 1
            WRITE(LOU,'(4E14.5)') ALFS, SGT, AFB
         ENDIF
1000  CONTINUE
      WRITE(LOU,*) 'JOIN SOLID'
C--
C  That's it.
C--
      STOP
      END
