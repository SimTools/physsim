      SUBROUTINE INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
#include "include/SMCONS.inc"
#include "include/SMPTAB.inc"
#include "include/SMCUPL.inc"
      REAL*8   ALF, ALFS, AMSW, AMSZ, AMSH, AMST
      REAL*8   BRT(12)
      EXTERNAL BLKDAT
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants.
C--
      xPI    = ACOS(-1.D0)
      x2PI   = 2*xPI
      x4PI   = 4*xPI
      xHPI   = xPI/2
      xSQRT2 = SQRT(2.0D0)
      xGV2PB = 3.8937966D8
      xGV2FB = 3.8937966D11
      xRD2DG = 180/xPI
      xDG2RD = xPI/180
C--
C  Initialize coupling constants.
C--
      xALF   = ALF
      xSIN2W = (1-AMSW/AMSZ)*(1+AMSW/AMSZ)
      xCOS2W = 1 - xSIN2W
      xALFS  = ALFS
C--
C  Initialize particle table.
C--
      AMW        = AMSW
      AMZ        = AMSZ
      AMH        = AMSH
      AMF(3,1,2) = AMST
C--
C  Print out modified parameters.
C--
      PRINT *, ' ****** INSMCN MODIFIED /SMCONS/ AND /SMPTAB/ ********'
      PRINT *, '             '
      PRINT *, '    ALF    = ', xALF
      PRINT *, '    ALFS   = ', xALFS
      PRINT *, '    S2W    = ', xSIN2W
      PRINT *, '    AMW    = ', AMW
      PRINT *, '    AMZ    = ', AMZ
      PRINT *, '    AMH    = ', AMH
      PRINT *, '    AMT    = ', AMF(3,1,2)
C--
C  Initialize W width.
C--
      COLOR  = 3
      QCD    = 1 + xALFS/xPI
      NMD    = 0
      GMWTOT = 0
C--Leptons.
      C   = 1
      VFF = 1
      DO 10 IG = 1, 3
         CALL GMW2FF(AMW,AMF(IG,1,1),AMF(IG,2,1),VFF,C,GM)
         GMWTOT   = GMWTOT + GM
         NMD      = NMD + 1
         BRT(NMD) = GMWTOT
10    CONTINUE
C--Quarks.
      C = COLOR*QCD
      DO 20 IGD = 1, 3
         DO 2 IGU = 1, 3
            VFF = VKM(IGU,IGD,2)
            CALL GMW2FF(AMW,AMF(IGU,1,2),AMF(IGD,2,2),VFF,C,GM)
            GMWTOT   = GMWTOT + GM
            NMD      = NMD + 1
            BRT(NMD) = GMWTOT
2        CONTINUE
20    CONTINUE
C--
C  Set branching fractions.
C--
      DO 25 IMD = 1, NMD
         BRVDK(IMD,1) = BRT(IMD)/GMWTOT
25    CONTINUE
C--
      PRINT *, '        '
      PRINT *, '    GMWTOT = ', GMWTOT
      PRINT *, '    BR(lp) = ', BRT(3)/GMWTOT
      PRINT *, '    BR(qk) = ', (BRT(12)-BRT(3))/GMWTOT
C--
C  Initialize Z width.
C--
      NMD    = 0
      GMZTOT = 0
C--Leptons.
      C = 1
      DO 30 IT = 1, 2
         T3F = 1.5D0 - IT
         QF  = 1.0D0 - IT
         DO 3 IG = 1, 3
            CALL GMZ2FF(AMZ,T3F,QF,C,AMF(IG,IT,1),GM)
            GMZTOT   = GMZTOT + GM
            NMD      = NMD + 1
            BRT(NMD) = GMZTOT
3        CONTINUE
30    CONTINUE
C--Quarks.
      C = COLOR*QCD
      DO 40 IT = 1, 2
         T3F = 1.5D0  - IT
         QF  = 5.D0/3 - IT
         DO 4 IG = 1, 3
            CALL GMZ2FF(AMZ,T3F,QF,C,AMF(IG,IT,2),GM)
            GMZTOT   = GMZTOT + GM
            NMD      = NMD + 1
            BRT(NMD) = GMZTOT
4        CONTINUE
40    CONTINUE
C--
C  Set branching fractions.
C--
      DO 45 IMD = 1, NMD
         BRVDK(IMD,2) = BRT(IMD)/GMZTOT
45    CONTINUE
C--
      PRINT *, '        '
      PRINT *, '    GMZTOT = ', GMZTOT
      PRINT *, '    BR(nu) = ', BRT(3)/GMZTOT
      PRINT *, '    BR(lp) = ', (BRT( 6)-BRT( 3))/GMZTOT
      PRINT *, '    BR(qu) = ', (BRT( 9)-BRT( 6))/GMZTOT
      PRINT *, '    BR(qd) = ', (BRT(12)-BRT( 9))/GMZTOT
C--
C  Initialize Higgs width.
C--
      AMSTU = AMF(3,2,1)
      AMSC  = AMF(2,1,2)
      AMSB  = AMF(3,2,2)
      CALL GMH2FF(AMSH,AMSTU,AMSZ,GMH2TU)
      CALL GMH2FF(AMSH,AMSC,AMSZ,GMH2C)
      CALL GMH2FF(AMSH,AMSB,AMSZ,GMH2B)
      CALL GMH2FF(AMSH,AMST,AMSZ,GMH2T)
      CALL GMH2VV(AMSH,AMSW,AMSZ,GMH2W)
      CALL GMH2VV(AMSH,AMSZ,AMSZ,GMH2Z)
      C      = COLOR
      GMH2C  = C*GMH2C
      GMH2B  = C*GMH2B
      GMH2T  = C*GMH2T
      GMH2W  = 2*GMH2W
      GMHTOT = GMH2TU + GMH2C + GMH2B + GMH2T + GMH2W + GMH2Z
C--
      PRINT *, '        '
      PRINT *, '    GMHTOT = ', GMHTOT
      IF ( GMHTOT.LE.0.D0 ) GMHTOT = 1.D-10
      PRINT *, '    BR(tu) = ', GMH2TU/GMHTOT
      PRINT *, '    BR(c)  = ', GMH2C/GMHTOT
      PRINT *, '    BR(b)  = ', GMH2B/GMHTOT
      PRINT *, '    BR(t)  = ', GMH2T/GMHTOT
      PRINT *, '    BR(W)  = ', GMH2W/GMHTOT
      PRINT *, '    BR(Z)  = ', GMH2Z/GMHTOT
C--
C  Initialize top width.
C--
      VFF = (VKM(3,3,2))**2
      CALL GMT2QW(AMST,AMSB,AMSW,VFF,GM)
      GMF(3,1,2) = GM
C--
      PRINT *, '        '
      PRINT *, '    GMTTOT = ', GMF(3,1,2)
      PRINT *, '        '
      PRINT *, ' *****************************************************'
C--
C  Initialize coupling constants.
C--
      CALL UHCUPL(ALF,xSIN2W, GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1,
     .            GW,GWWZ,GWWA)
C--
C  That's it.
C--
      RETURN
      END

