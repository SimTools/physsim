C
C=======================================================================
C  Subroutines necessary to run this job.
C=======================================================================
C
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
C  *INCLUDE SGAFB0
C  *INCLUDE SGAFBF
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
C
CC**********************************************************************
C
C*
C*  92/10/23  K.Fujii        This program plots AFB and SIG for
C*                           various VTB values.
C*
      PROGRAM pltsgav

      IMPLICIT REAL*8  ( A-H, O-Z )
 
      CHARACTER*8   JOIN(0:3)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOU      /20/
C--
C  Allocate output data set.
C--
C      CALL UALCPS(LOU,'TKSF.@@.@SGSPA@.VTB2.IR.TDR','RENEW','FB',30,IRT)
C     CALL UALCPS(LOU,'TKSF.@@.@SGSPA@.VTB2.BM.TDR','RENEW','FB',30,IRT)
C--
C  Set constants.
C--
      ALP  = 1/128.D0
      SN2W = 0.23D0
      AMSZ = 91.17D0
      AMSW = 80.0D0
      AMSB =  5.0D0
C--
C  Set top parameters.
C--
      ALPS   = 0.12D0
C>>>
      AMST   = 150.D0
C>>>
      VTB2MN = 0.8D0
      VTB2MX = 1.2D0
      NVTB2  = 2
      DVTB2  = (VTB2MX-VTB2MN)/NVTB2
C>>>
      NVTB2  = 0
      VTB2MN = 1
C>>>
C--
C  Initialize radiation function.
C--
      CALL ISRINT(0,2*AMST)
C--
C  Energy range.
C--
      EMN  = -8
      EMX  = +4
      NE   = 60
C>>>
C     EMN  = -4
C     EMX  = +1
C     NE   = 200
C>>>
      DE   = (EMX-EMN)/NE
C--
C  Prepare TOPDRAW files.
C--
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''( RQQ'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''NEWFRAME'')')
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''SET INTENSITY 4'')')
      WRITE(LOU,'(''SET WINDOW X  2. 12. Y 1.8  9.2'')')
      WRITE(LOU,'(''SET LIMITS X '',2F9.0,'' Y 0 3.0'')') EMN, EMX
      WRITE(LOU,'(''SET TITLE SIZE 3'')')
      WRITE(LOU,'(''TITLE 0.5 5.0 ANGLE 90 SIZE 4 ''''S0tt1(pb)'')')
      WRITE(LOU,'(''CASE                          ''''GX  X    '')')
      WRITE(LOU,'(''TITLE 6.0 0.8 SIZE 4 ''''E (GeV)'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOU,'(''( AMST = '',F8.3)') AMST
      WRITE(LOU,'(''( ALPS = '',F8.3)') ALPS
      WRITE(LOU,'(''SET ORDER X Y DUMMY'')')
C--
C  Loop over energies.
C--
      MODE = 0
      DO 2000 IV = 0, NVTB2
         VTB2 = VTB2MN + IV*DVTB2
         WRITE(LOU,'(''( VTB2 = '',F8.3)') VTB2
         DO 200 IE = 0, NE
            E = EMN + IE*DE
            CALL SGAFBF(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,
     .                  E,SIGTOT,DLTAFB)
            WRITE(LOU,*) E, SIGTOT, DLTAFB
            PRINT *,    E, SIGTOT, DLTAFB
            MODE = 3
200      CONTINUE
         WRITE(LOU,'(''JOIN '',A)') JOIN(IV)
         MODE = 2
2000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
