      SUBROUTINE INSSCN(ALF,ALFS,WM,ZM,AM0,AMU,AM2,TANB,AMA)
 
      IMPLICIT    REAL*4 ( A-H, O-Z )
#include "../ss/include/SSCONS.inc"
#include "../ss/include/SSPTAB.inc"
#include "../ss/include/SSCUPL.inc"
      REAL*4   ALF, ALFS, WM, ZM, AM0, AMU, AM2, TANB, AMA
      REAL*4   FM(3)
C
C========< Entry Point >================================================
C
C--
C  Initialize SUSY parameters.
C--
      xAM0   = AM0
      xAMU   = AMU
      xAM2   = AM2
      xTANB  = TANB
      xAMA   = AMA
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
C--
      GMSH(1) = 0.5
      GMSH(2) = 0.5
      GMSH(3) = 0.5
      GMSH(4) = 0.5
C--
C  Print out modified parameters.
C--
      PRINT *, '             '
      PRINT *, ' ****** INSSCN MODIFIED /SSCONS/ AND /SSPTAB/ ********'
      PRINT *, '             '
      PRINT *, '    AM0    = ', xAM0
      PRINT *, '    AMU    = ', xAMU
      PRINT *, '    AM2    = ', xAM2
      PRINT *, '    TANB   = ', xTANB
      PRINT *, '    AMA    = ', xAMA
      PRINT *, '        '
      PRINT *, '    AMSNL  = ', SFM(1)
      PRINT *, '    AMSEL  = ', SFM(2)
      PRINT *, '    AMSER  = ', SFM(3)
      PRINT *, '    AMSUL  = ', SFM(4)
      PRINT *, '    AMSUR  = ', SFM(5)
      PRINT *, '    AMSDL  = ', SFM(6)
      PRINT *, '    AMSDR  = ', SFM(7)
      PRINT *, '        '
      PRINT *, '    AMSW1  = ', SWM(1)
      PRINT *, '    AMSW2  = ', SWM(2)
      PRINT *, '        '
      PRINT *, '    AMSZ1  = ', SZM(1)
      PRINT *, '    AMSZ2  = ', SZM(2)
      PRINT *, '    AMSZ3  = ', SZM(3)
      PRINT *, '    AMSZ4  = ', SZM(4)
      PRINT *, '        '
      PRINT *, '    GMSNL  = ', GMSF(1)
      PRINT *, '    GMSEL  = ', GMSF(2)
      PRINT *, '    GMSER  = ', GMSF(3)
      PRINT *, '    GMSUL  = ', GMSF(4)
      PRINT *, '    GMSUR  = ', GMSF(5)
      PRINT *, '    GMSDL  = ', GMSF(6)
      PRINT *, '    GMSDR  = ', GMSF(7)
      PRINT *, '        '
      PRINT *, '    GMSW1  = ', GMSW(1)
      PRINT *, '    GMSW2  = ', GMSW(2)
      PRINT *, '        '
      PRINT *, '    GMSZ1  = ', GMSZ(1)
      PRINT *, '    GMSZ2  = ', GMSZ(2)
      PRINT *, '    GMSZ3  = ', GMSZ(3)
      PRINT *, '    GMSZ4  = ', GMSZ(4)
      PRINT *, '        '
      PRINT *, ' *****************************************************'
C--
C  That's it.
C--
      RETURN
      END
