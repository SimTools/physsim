C  C
C  C====================================================================
C  C  Subroutines necessary to run this job.
C  C====================================================================
C  C
C  C--
C  C  QCD corrections.
C  C--
C  *INCLUDE GTCEFF
C  *INCLUDE GMTWMB
C  *INCLUDE QCDPOT
C  *INCLUDE ADJ1
C  *INCLUDE ADJ2
C  *INCLUDE QCDS1
C  *INCLUDE QCDS2
C  *INCLUDE QCDS3
C  *INCLUDE BNDCND
C  *INCLUDE RKUTAS
C  *INCLUDE QCDP1
C  *INCLUDE QCDP2
C  *INCLUDE QCDP3
C  *INCLUDE RKUTAP
C  *INCLUDE DELFN
C  *INCLUDE RUNWID
C  *INCLUDE LNT0FI
C  *INCLUDE LNT1FI
C  *INCLUDE NTGRD1
C  *INCLUDE NTGRD2
C  *INCLUDE NTGRD3
C  *INCLUDE LNT0WA
C  *INCLUDE PHS0
C  *INCLUDE LNT2WA
C  *INCLUDE PHS2
C  *INCLUDE TBWCORR
C  *INCLUDE FX
C  *INCLUDE GTQCDP
C  *INCLUDE GAMMA0
C  C--
C  C  Cross section and asymmetry.
C  C--
C  *INCLUDE SGAFBF
C  *INCLUDE SGAFB0
C  *INCLUDE GRQSGA
C  *INCLUDE GRQINT
C  *INCLUDE GRQADJ
C  *INCLUDE GRQSUM
C  C--
C  C  ISR and beam effects.
C  C       RSHDIS0 : ISR only.
C  C       RSHDIS  : ISR + Beam effects.
C  C--
C  *INCLUDE UDSRCH
C  C*INCLUDE RSHDIS
C  *INCLUDE RSHDIS0
C  C--
C  C  E_peak hunting.
C  C--
C  *INCLUDE EPKHNT
C  *INCLUDE SGTTHR2
C  *INCLUDE RKUTTA2
C  *INCLUDE EI
C  *INCLUDE F1
C  *INCLUDE FH
C  *INCLUDE GTGAMT
C  *INCLUDE SETPRMC
CC**********************************************************************
C*
C* (Purpose)
C*     This program calculates the FB-asymmetry as a function of
C*     DE = E - E_1S for various alpha_s values.
C* (Update Record)
C*    93/04/19  K.Fujii                Original version w/ ISR.
C*
CC**********************************************************************

      PROGRAM pltfb1sea
 
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
C     CALL UALCPS(LOU,'TKSF.@@.@FBDE@.ALFS.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@FBDE@.ALFS.BM.TDR','RENEW','FB',30,IRT)
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
      NDE    =  25
      DEMN   = -1.D0
      DEMX   = +4.D0
      DDE    = (DEMX-DEMN)/NDE
C--
C  Initialize ISR.
C--
      CALL ISRINT(0,2*AMT)
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
               CALL SGAFBF(MODE2,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
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
