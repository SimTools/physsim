CINCLUDE INSUSY
CINCLUDE COUPLSUB
CINCLUDE DCY2BODY
CINCLUDE GTSFMS
CINCLUDE SF2BR
CINCLUDE INOMIX
CINCLUDE HGSMIX
CINCLUDE SWPELM
 
      IMPLICIT    REAL*4 ( A-H, O-Z )
C=EXPAND 'TKSF.G#XCXC.FORT(@SSCONS)'
C=EXPAND 'TKSF.G#XCXC.FORT(@SSPTAB)'
      COMMON /SSPTAB/ SWM, SZM, SFM, SHM, GMSW, GMSZ, GMSF, GMSH
      REAL*4          SWM(2), SZM(4), SFM(7), SHM(4),
     .                GMSW(2), GMSZ(4), GMSF(7), GMSH(4)
C=EXPAND 'TKSF.G#XCXC.FORT(@SSCUPL)'
      COMMON /SSCUPL/ GNCW, GCCZ, GNNZ,
     .                GCESNL, GCNSEL, GCNSER,
     .                GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     .                GNNSNL, GNESEL, GNESER,
     .                GNUSUL, GNUSUR, GNDSDL, GNDSDR
C--
      COMPLEX*8  GNCW(2,4,2), GCCZ(2,2,2), GNNZ(2,4,4)
      COMPLEX*8  GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     .           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*8  GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     .           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
C--
      REAL*4    ALF, ALFS, WM, ZM, AM0, AMU, AM2, TANB, AMA
      REAL*4    FM(3)
      COMPLEX*8 GC(2)
C
C========< Entry Point >================================================
C
C--
C  Initialize SUSY parameters.
C--
      AM0   =  20
      AMU   = 400
      AM2   = 250
      TANB  = 2
      AMA   = 600
      WM    = 80
      ZM    = 91.17
      ALF   = 1./128.
      ALFS  = 0.12
C--
C  Initialize particle table and coupling constants.
C  This version ignors Yukawa couplings.
C--
      S2W   = ( 1 - WM/ZM )*( 1 + WM/ZM )
      FM(1) = 0
      FM(2) = 0
      FM(3) = 0
C--
      CALL INSUSY(AM0,AMU,AM2,TANB,AMA,ALF,ALFS,S2W,ZM,FM,
     .            SFM,SWM,SZM,GMSF,GMSW,GMSZ,
     .            GNCW,GCCZ,GNNZ,
     .            GCESNL,GCNSEL,GCNSER,
     .            GCUSDL,GCUSDR,GCDSUL,GCDSUR,
     .            GNNSNL,GNESEL,GNESER,
     .            GNUSUL,GNUSUR,GNDSDL,GNDSDR )
C>>>
      WRITE(20,*) '             '
      WRITE(20,*) ' ****** Coupling Constants in /SSCUPL/ ************'
      WRITE(20,*) '             '
      WRITE(20,*) '  GNCW((*,1,1) = ', GNCW(1,1,1), GNCW(2,1,1)
      WRITE(20,*) '             '
      WRITE(20,*) '  GCESNL(*,1) = ', GCESNL(1,1), GCESNL(2,1)
      WRITE(20,*) '  GCNSEL(*,1) = ', GCNSEL(1,1), GCNSEL(2,1)
      WRITE(20,*) '  GCNSER(*,1) = ', GCNSER(1,1), GCNSER(2,1)
      WRITE(20,*) '  GCUSDL(*,1) = ', GCUSDL(1,1), GCUSDL(2,1)
      WRITE(20,*) '  GCUSDR(*,1) = ', GCUSDR(1,1), GCUSDR(2,1)
      WRITE(20,*) '  GCDSUL(*,1) = ', GCDSUL(1,1), GCDSUL(2,1)
      WRITE(20,*) '  GCDSUR(*,1) = ', GCDSUR(1,1), GCDSUR(2,1)
      WRITE(20,*) ' '
      WRITE(20,*) '  GNNSNL(*,1) = ', GNNSNL(1,1), GNNSNL(2,1)
      WRITE(20,*) '  GNESEL(*,1) = ', GNESEL(1,1), GNESEL(2,1)
      WRITE(20,*) '  GNESER(*,1) = ', GNESER(1,1), GNESER(2,1)
      WRITE(20,*) '  GNUSUL(*,1) = ', GNUSUL(1,1), GNUSUL(2,1)
      WRITE(20,*) '  GNUSUR(*,1) = ', GNUSUR(1,1), GNUSUR(2,1)
      WRITE(20,*) '  GNDSDL(*,1) = ', GNDSDL(1,1), GNDSDL(2,1)
      WRITE(20,*) '  GNDSDR(*,1) = ', GNDSDR(1,1), GNDSDR(2,1)
C--
C  Print out modified parameters.
C--
      WRITE(20,*) '             '
      WRITE(20,*) ' ****** INSSCN MODIFIED /SSCONS/ AND /SSPTAB/ *****'
      WRITE(20,*) '             '
      WRITE(20,*) '    AM0    = ', AM0
      WRITE(20,*) '    AMU    = ', AMU
      WRITE(20,*) '    AM2    = ', AM2
      WRITE(20,*) '    TANB   = ', TANB
      WRITE(20,*) '    AMA    = ', AMA
      WRITE(20,*) '        '
      WRITE(20,*) '    AMSNL  = ', SFM(1)
      WRITE(20,*) '    AMSEL  = ', SFM(2)
      WRITE(20,*) '    AMSER  = ', SFM(3)
      WRITE(20,*) '    AMSUL  = ', SFM(4)
      WRITE(20,*) '    AMSUR  = ', SFM(5)
      WRITE(20,*) '    AMSDL  = ', SFM(6)
      WRITE(20,*) '    AMSDR  = ', SFM(7)
      WRITE(20,*) '        '
      WRITE(20,*) '    AMSW1  = ', SWM(1)
      WRITE(20,*) '    AMSW2  = ', SWM(2)
      WRITE(20,*) '        '
      WRITE(20,*) '    AMSZ1  = ', SZM(1)
      WRITE(20,*) '    AMSZ2  = ', SZM(2)
      WRITE(20,*) '    AMSZ3  = ', SZM(3)
      WRITE(20,*) '    AMSZ4  = ', SZM(4)
      WRITE(20,*) '        '
      WRITE(20,*) '    GMSNL  = ', GMSF(1)
      WRITE(20,*) '    GMSEL  = ', GMSF(2)
      WRITE(20,*) '    GMSER  = ', GMSF(3)
      WRITE(20,*) '    GMSUL  = ', GMSF(4)
      WRITE(20,*) '    GMSUR  = ', GMSF(5)
      WRITE(20,*) '    GMSDL  = ', GMSF(6)
      WRITE(20,*) '    GMSDR  = ', GMSF(7)
      WRITE(20,*) '        '
      WRITE(20,*) '    GMSW1  = ', GMSW(1)
      WRITE(20,*) '    GMSW2  = ', GMSW(2)
      WRITE(20,*) '        '
      WRITE(20,*) '    GMSZ1  = ', GMSZ(1)
      WRITE(20,*) '    GMSZ2  = ', GMSZ(2)
      WRITE(20,*) '    GMSZ3  = ', GMSZ(3)
      WRITE(20,*) '    GMSZ4  = ', GMSZ(4)
      WRITE(20,*) '        '
      PRINT *, ' **************************************************'
C--
      PRINT *, '*** Chargino Decay Widths ***'
      GAMTOT = 0
      CALL GAMFFV(SWM(1),SZM(1),WM,GNCW(1,1,1),GAM)
      PRINT *, ' GAMMA_X+ = ', GAM
      GAMTOT = GAMTOT + GAM
C--
      GAMSNL = 0
      GC(1) = CONJG(GCESNL(2,1))
      GC(2) = CONJG(GCESNL(1,1))
      AME   = 0.511E-3
      CALL GAMFFS(SWM(1),AME,SFM(1),GC,GAM)
      PRINT *, ' GAMMA_SNE = ', GAM
C--
      AMMU   = 0.106
      CALL GAMFFS(SWM(1),AMMU,SFM(1),GC,GAM)
      PRINT *, ' GAMMA_SNM = ', GAM
      GAMSNL = GAMSNL + GAM
C--
      AMTU   = 1.777
      CALL GAMFFS(SWM(1),AMTU,SFM(1),GC,GAM)
      PRINT *, ' GAMMA_SNT = ', GAM
      GAMSNL = GAMSNL + GAM
      PRINT *, ' GAMMA_SNL = ', GAMSNL
      GAMTOT = GAMTOT + GAMSNL
C--
      GC(1) = GCNSEL(1,1)
      GC(2) = GCNSEL(2,1)
      CALL GAMFFS(SWM(1),0.,SFM(2),GC,GAM)
      PRINT *, ' GAMMA_SLN = ', GAM
      PRINT *, ' GAMMA_SLN = ', 3*GAM
      GAMTOT = GAMTOT + 3*GAM
      PRINT *, ' GAMMA_TOT = ',  GAMTOT
C--
      PRINT *, '*** Slepton Decay Widths ***'
      GC(1) = GNESEL(1,1)
      GC(2) = GNESEL(2,1)
      CALL GAMSFF(SFM(2),SZM(1),AME,GC,GAM)
      PRINT *, ' GAMMA_SEL = ',  GAM
C--
      CALL GAMSFF(SFM(2),SZM(1),AMMU,GC,GAM)
      PRINT *, ' GAMMA_SML = ',  GAM
C--
      CALL GAMSFF(SFM(2),SZM(1),AMTU,GC,GAM)
      PRINT *, ' GAMMA_STL = ',  GAM
C--
      GC(1) = GNNSNL(1,1)
      GC(2) = GNNSNL(2,1)
      CALL GAMSFF(SFM(1),0.,SZM(1),GC,GAM)
      PRINT *, ' GAMMA_SNL = ',  GAM
      CALL SCF1F2(SFM(1),0.,SZM(1),GC(1),GC(2),GAM)
      PRINT *, ' GAMMA_SNL = ',  GAM
C--
      PRINT *, '*** Squark Decay Widths ***'
      GC(1) = GNDSDL(1,1)
      GC(2) = GNDSDL(2,1)
      CALL GAMSFF(SFM(2),SZM(1),0.,GC,GAM)
      PRINT *, ' GAMMA_SDL = ',  GAM
C--
      GC(1) = GNUSUL(1,1)
      GC(2) = GNUSUL(2,1)
      CALL GAMSFF(SFM(2),0.,SZM(1),GC,GAM)
      PRINT *, ' GAMMA_SUL = ',  GAM
C--
C  That's it.
C--
      STOP
      END
