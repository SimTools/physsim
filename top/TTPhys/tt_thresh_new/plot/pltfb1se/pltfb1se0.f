C
C=======================================================================
C  Subroutines necessary to run this job.
C=======================================================================
C
C--
C  QCD corrections.
C--
CINCLUDE GTCEFF
CINCLUDE GMTWMB
CINCLUDE QCDPOT
CINCLUDE ADJ1
CINCLUDE ADJ2
CINCLUDE QCDS1
CINCLUDE QCDS2
CINCLUDE QCDS3
CINCLUDE BNDCND
CINCLUDE RKUTAS
CINCLUDE QCDP1
CINCLUDE QCDP2
CINCLUDE QCDP3
CINCLUDE RKUTAP
CINCLUDE DELFN
CINCLUDE RUNWID
CINCLUDE LNT0FI
CINCLUDE LNT1FI
CINCLUDE NTGRD1
CINCLUDE NTGRD2
CINCLUDE NTGRD3
CINCLUDE LNT0WA
CINCLUDE PHS0
CINCLUDE LNT2WA
CINCLUDE PHS2
CINCLUDE TBWCORR
CINCLUDE FX
CINCLUDE GTQCDP
CINCLUDE GAMMA0
C--
C  Cross section and asymmetry.
C--
CINCLUDE SGAFB0
CINCLUDE GRQSGA
CINCLUDE GRQINT
CINCLUDE GRQADJ
CINCLUDE GRQSUM
C--
C  E_peak hunting.
C--
CINCLUDE EPKHNT
CINCLUDE SGTTHR2
CINCLUDE RKUTTA2
CINCLUDE EI
CINCLUDE F1
CINCLUDE FH
CINCLUDE GTGAMT
CINCLUDE SETPRMC
CC**********************************************************************
C*
C* (Purpose)
C*     This program calculates the FB-asymmetry as a function of
C*     DE = E - E_1S for various alpha_s values.
C* (Update Record)
C*    93/04/19  K.Fujii                Original version w/ ISR.
C*
CC**********************************************************************

      PROGRAM pltfb1se0
 
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
C     CALL UALCPS(LOU,'TKSF.@@.@FBDE@.ALFS.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@FBDE@.ALFS.NI.TDR','RENEW','FB',30,IRT)
C--
C  Initialize constants.
C--
      VTB2   = 1.0D0
      AMH    = 1.D10
      BTH    = 0.D10
C--
      ALFSMN = 0.11D0
      ALFSMX = 0.13D0
      NA     = 2
      DA     = (ALFSMX-ALFSMN)/NA
C
      ALFSMN = 0.12D0
      NA     = 0
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
C  Energy range.
C--
      NDE    =  10
      DEMN   = -1.D0
      DEMX   = +4.D0
      DDE    = (DEMX-DEMN)/NDE
C--
C  Topdraw commands.
C--
      WRITE(LOU,*) '( '
      WRITE(LOU,*) '( m_t      = ', REAL(AMT)
      WRITE(LOU,*) '( !V_tb!^2 = ', REAL(VTB2)
      WRITE(LOU,*) '( No Higgs '
      WRITE(LOU,*) '( '
      WRITE(LOU,*) 'SET LIMITS X ', DEMN, DEMX, ' Y -0.05 0.15'
      WRITE(LOU,*) 'SET ORDER X DUMMY DUMMY Y'
C--
C  Loop over alpha_s.
C--
      MODE1 = 0
      MODE2 = 0
      DO 1000 IA = 0, NA
         ALFS = ALFSMN + IA*DA
         WRITE(LOU,*) '( alpha_s  = ', REAL(ALFS)
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
               CALL SGAFB0(MODE2,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                     E,SGT,AFB)
               WRITE(LOU,'(4E14.5)') DE, E, SGT, AFB
               MODE2 = 3
10          CONTINUE
            WRITE(LOU,*) 'JOIN ', JOIN(IA)
         ENDIF
         MODE1 = 0
         MODE2 = 1
1000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
