C*
C*  92/10/23  K.Fujii        This program plots AFB and SIG for
C*                           various alpha_s values.
C*
      PROGRAM pltsgaa 
 
      IMPLICIT REAL*8  ( A-H, O-Z )
 
      CHARACTER*8   JOIN(0:3)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOU      /20/
C--
C  Allocate output data set.
C--
C      CALL UALCPS(LOU,'TKSF.@@.@SGSPA@.ALFS.TDR','RENEW','FB',30,IRT)
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
      AMST = 150.D0
      VTB2 = 1.D0
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
            CALL SGAFBF(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,
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
