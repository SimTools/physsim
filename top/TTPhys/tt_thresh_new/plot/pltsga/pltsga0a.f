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
C
C=======================================================================
C  Main programs.
C=======================================================================
C
C*
C*  92/10/23  K.Fujii        This program plots AFB and SIG for
C*                           various alpha_s values.
C*
 
      IMPLICIT REAL*8  ( A-H, O-Z )
 
      CHARACTER*8   JOIN(0:3)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOU      /21/
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
      ASMN = 0.11
      ASMX = 0.13
      NA   = 2
      DA   = (ASMX-ASMN)/NA
C>>>
      ASMN = 0.12
      NA   = 0
C>>>
C--
C     AMST = 150.D0
      AMST = 170.D0
C>>>
      VTB2 = 1
C     VTB2 = 0.1
C>>>
C--
C  Energy range.
C--
C>>>
      EMN  = -8
      EMX  = +4
      NE   = 60
C>>>
C      EMN  = -4
C      EMX  = +1
C      NE   = 200
C>>>
C     EMN  = -0.6
C     EMX  = +1
C     NE   = 16
C>>>
      DE   = (EMX-EMN)/NE
C--
C  Allocate output data set.
C--
C     CALL UALCPS(LOU,'TKSF.@@.@SGFB@.VT01.M5S','RENEW','FB',30,IRT)
C     CALL UALCPS(LOU,'TKSF.@@.@SGFB@.VT01.NOM','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@SGFB@.VT01.P5S','RENEW','FB',30,IRT)
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
      WRITE(LOU,'(''( VTB2 = '',F8.3)') VTB2
      WRITE(LOU,'(''SET ORDER X Y DUMMY'')')
C--
C  Loop over energies.
C--
      DO 2000 IA = 0, NA
         ALPS = ASMN + IA*DA
         WRITE(LOU,'(''( ALPS = '',F8.3)') ALPS
         MODE = 0
         DO 200 IE = 0, NE
            E = EMN + IE*DE
            CALL SGAFB0(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,
     .                  E,SIGTOT,DLTAFB)
            WRITE(LOU,*) E, SIGTOT, DLTAFB
C           PRINT *,    E, SIGTOT, DLTAFB
            MODE = 3
200      CONTINUE
         WRITE(LOU,'(''JOIN '',A)') JOIN(IA)
2000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
