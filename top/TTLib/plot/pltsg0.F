CINCLUDE EI
C>>>
C*INCLUDE POTRCH
CINCLUDE POTQCD
C*INCLUDE POTCPL
C*INCLUDE POTCLM
C*INCLUDE POTMTN
C*INCLUDE POTHIG
C>>>
CINCLUDE RKUTTA3
CINCLUDE SETPRM
CINCLUDE GTGAMT
 
C
C  sigma_tt0 for various alpha_s
C
C
      IMPLICIT REAL*8 ( A-H, O-Z )
      PARAMETER     ( NPT = 280 )
      REAL*8        RDATA(2,0:NPT,0:10)
      CHARACTER*8   JOIN(0:3)
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOUNIT   / 21 /
C--
C  Initialize parameters.
C--
      PI   = ACOS(-1.D0)
      ALP  = 1./128.
      SN2W = 0.23
      AMSZ = 91.18
      GAMZ = 2.5
      AMSW = 80.0
      AMSB = 5.0
      AMST = 175
      VFF  = 1.0
      AMSH = 1.D8
      BETH = 0
C--
C  Alpha_s range.
C--
      NA     = 2
      ALFSMN = 0.11
      ALFSMX = 0.13
      DA     = (ALFSMX-ALFSMN)/NA
C>>>
      NA     = 0
      ALFSMN = 0.12
C     NA     = 0
C     ALFSMN = 0.16
C--
C  Energy range.
C--
      RSMN = 2*AMST - 20
      RSMX = 2*AMST + 50
      DRS  = (RSMX-RSMN)/NPT
C--
C  Production cross section.
C--
      DO 2000 IA = 0, NA
         ALFS = ALFSMN + IA*DA
         MODE = 0
         DO 200 I = 0, NPT
            RS = RSMN + DRS*I
            CALL SGTT0(MODE,ALFS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                 AMST, VFF, AMSH,BETH,   RS,SG)
            RDATA(1,I,IA) = RS
            RDATA(2,I,IA) = SG
C>>>
C        PRINT *, ' RS, SG = ', RS, SG
C>>>
            MODE       = 2
200      CONTINUE
2000  CONTINUE
C--
C  Prepare TOPDRAW files.
C--
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''( RQQ'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''NEWFRAME'')')
      WRITE(LOUNIT,'(''SET INTENSITY 4'')')
      WRITE(LOUNIT,'(''SET FONT DUPLEX'')')
      WRITE(LOUNIT,'(''SET WINDOW X  2. 12. Y 1.8  9.2'')')
      WRITE(LOUNIT,'(''SET LIMITS X '',2F9.0,'' Y 0 1500'')') RSMN, RSMX
      WRITE(LOUNIT,'(''SET TITLE SIZE 3'')')
      WRITE(LOUNIT,'(''TITLE 0.5 5.0 ANGLE 90 SIZE 4 ''''S0tt1(pb)'')')
      WRITE(LOUNIT,'(''CASE                          ''''GX  X    '')')
      WRITE(LOUNIT,'(''TITLE 6.0 0.8 SIZE 4 ''''2s0O (GeV)'')')
      WRITE(LOUNIT,'(''CASE                       ''''M UD'')')
      WRITE(LOUNIT,'(''('')')
      WRITE(LOUNIT,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOUNIT,'(''SET ORDER X Y'')')
C--
      DO 3000 IA = 0, NA
         DO 300 I = 0, NPT
            WRITE(LOUNIT,'(5X,5F12.5)') RDATA(1,I,IA), RDATA(2,I,IA)
300      CONTINUE
         WRITE(LOUNIT,'(''JOIN '',A8)') JOIN(IA)
3000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
