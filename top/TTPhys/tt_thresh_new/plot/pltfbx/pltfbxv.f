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
C*     !V_tb!^2 at a given delta_E value.
C* (Update Record)
C*    93/04/19  K.Fujii                Original version w/ ISR.
C*
CC**********************************************************************

      PROGRAM pltfbxv
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      DATA LOU /21/
C
C========< Entry Point >================================================
C
C--
C  Allocate output file.
C--
C     CALL UALCPS(LOU,'TKSF.@@.@FBDX@.VTB2.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@FBDX@.VTB2.BM.TDR','RENEW','FB',30,IRT)
C--
C  Initialize constants.
C--
      ALFS   = 0.12D0
      AMH    = 1.D10
      BTH    = 0.D10
C--
C>>>
      VTB2MN = 0.60D0
      VTB2MX = 1.40D0
      NV     = 1
      DV     = (VTB2MX-VTB2MN)/NV
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
      WRITE(LOU,*) '( alpha_s  = ', REAL(ALFS)
      WRITE(LOU,*) '( No Higgs '
      WRITE(LOU,*) '( '
      WRITE(LOU,*) 'SET LIMITS X ', VTB2MN, VTB2MX, ' Y -0.05 0.15'
      WRITE(LOU,*) 'SET ORDER X DUMMY Y'
C--
C  Loop over alpha_s.
C--
      MODE2 = 0
      DO 1000 IV = 0, NV
         VTB2 = VTB2MN + IV*DV
         WRITE(LOU,*) '( !V_tb!^2 = ', REAL(VTB2)
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
            WRITE(LOU,'(4E14.5)') VTB2, SGT, AFB
         ENDIF
1000  CONTINUE
      WRITE(LOU,*) 'JOIN SOLID'
C--
C  That's it.
C--
      STOP
      END
