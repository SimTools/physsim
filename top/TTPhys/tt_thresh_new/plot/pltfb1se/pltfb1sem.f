CC**********************************************************************
C*
C* (Purpose)
C*     This program calculates the FB-asymmetry as a function of
C*     DE = E - E_1S for various AMT values.
C* (Update Record)
C*    93/04/19  K.Fujii                Original version w/ ISR.
C*
CC**********************************************************************

      PROGRAM pltfb1sem
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      CHARACTER*8   JOIN(0:3)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOU /21/
C
C========< Entry Point >================================================
C
C--
C  Allocate output file.
C--
C     CALL UALCPS(LOU,'TKSF.@@.@FBDE@.AMST.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@FBDE@.AMST.BM.TDR','RENEW','FB',30,IRT)
C--
C  Initialize constants.
C--
      ALFS   = 0.12D0
      VTB2   = 1
      AMH    = 1.D10
      BTH    = 0.D10
C--
C     AMTMN  = 149.D0
C     AMTMX  = 151.D0
      AMTMN  = 169.D0
      AMTMX  = 171.D0
      NM     = 2
      DM     = (AMTMX-AMTMN)/NM
C--
      EPK    = -2.5D0
C--
      ALF    = 1/128.D0
      S2W    = 0.23D0
      AMZ    = 91.17D0
      AMW    = 80.0D0
      AMB    =  5.0D0
C--
C  Energy range.
C--
      NDE    =  25
      DEMN   = -1.D0
      DEMX   = +4.D0
      DDE    = (DEMX-DEMN)/NDE
C--
C  Initialize ISR.
C--
      CALL ISRINT(0,AMTMN+AMTMX)
C--
C  Topdraw commands.
C--
      WRITE(LOU,*) '( '
      WRITE(LOU,*) '( !V_tb!^2 = ', REAL(VTB2)
      WRITE(LOU,*) '( alpha_s  = ', REAL(ALFS)
      WRITE(LOU,*) '( No Higgs '
      WRITE(LOU,*) '( '
      WRITE(LOU,*) 'SET LIMITS X ', DEMN, DEMX, ' Y -0.05 0.15'
      WRITE(LOU,*) 'SET ORDER X DUMMY DUMMY Y'
C--
C  Loop over m_t.
C--
      MODE1 = 0
      MODE2 = 0
      DO 1000 IM = 0, NM
         AMT  = AMTMN + IM*DM
         WRITE(LOU,*) '( m_t      = ', REAL(AMT)
         CALL EPKHNT( MODE1, ALF,S2W,AMZ,AMW,AMB,
     .                       AMT,ALFS,VTB2,AMH,BTH, EPK,SGPK, NTRY )
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
            DO 10 IDE = 0, NDE
               DE = DEMN + IDE*DDE
               E     = EPK + DE
               CALL SGAFBF(MODE2,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                     E,SGT,AFB)
               WRITE(LOU,'(4E14.5)') DE, E, SGT, AFB
               MODE2 = 3
10          CONTINUE
            WRITE(LOU,*) 'JOIN ', JOIN(IM)
         ENDIF
         MODE1 = 1
         MODE2 = 2
1000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
