      SUBROUTINE FULCXD(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SSCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SMCUPL.inc"
#include "SSCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2(0:3)
C--
      COMPLEX*8  SPINOR(6,20), WVCT(6), SCAL(3), AMP(0:3)
      COMPLEX*8  GWX(2,2), GCFSF(2,2,2), GNFSF(2,2,2)
      REAL   *4  AMSSF(2,2), GAMSF(2,2)
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL     = 1
         AMSW      = AMW
         GAMW      = GMWTOT
C>>>
C        QCD       =  1
         QCD       =  1 + xALFS/xPI
C>>>
         GWX(1,1)  = CONJG(GNCW(1,1,1))
         GWX(2,1)  = CONJG(GNCW(2,1,1))
         GWX(1,2)  = GNCW(1,1,1)
         GWX(2,2)  = GNCW(2,1,1)
C>>>
CTBW          GWX(1,1)  = GWF(1)
CTBW          GWX(2,1)  = GWF(2)
C>>>
      ENDIF
C--
C  Set masses and couplings depending on final states.
C--
      CALL RDOTxx(PV(0,1),PV(0,1),AMSX)
      AMSX  = SQRT(AMSX)
      CALL RDOTxx(PV(0,2),PV(0,2),AMSX0)
      AMSX0 = SQRT(AMSX0)
C--
      AM1   = AMSX
      AM2   = AMSX0
      AM3   = AMF(IDP(1,3),IDP(2,3),IDP(3,3))
      AM4   = AMF(IDP(1,4),IDP(2,4),IDP(3,4))
C--
      IF ( IDP(3,3).EQ.1 ) THEN
         AMSSF(  1,1) = SFM(1)
         GAMSF(  1,1) = GMSF(1)
         GCFSF(1,1,1) = GCESNL(1,1)
         GCFSF(2,1,1) = GCESNL(2,1)
         GNFSF(1,1,1) = CONJG(GNNSNL(2,1))
         GNFSF(2,1,1) = CONJG(GNNSNL(1,1))
         AMSSF(  2,1) = SFM(2)
         GAMSF(  2,1) = GMSF(2)
         GCFSF(1,2,1) = CONJG(GCNSEL(2,1))
         GCFSF(2,2,1) = CONJG(GCNSEL(1,1))
         GNFSF(1,2,1) = GNESEL(1,1)
         GNFSF(2,2,1) = GNESEL(2,1)
      ELSE
         AMSSF(  1,1) = SFM(4)
         GAMSF(  1,1) = GMSF(4)
         GCFSF(1,1,1) = GCDSUL(1,1)
         GCFSF(2,1,1) = GCDSUL(2,1)
         GNFSF(1,1,1) = CONJG(GNUSUL(2,1))
         GNFSF(2,1,1) = CONJG(GNUSUL(1,1))
         AMSSF(  2,1) = SFM(6)
         GAMSF(  2,1) = GMSF(6)
         GCFSF(1,2,1) = CONJG(GCUSDL(2,1))
         GCFSF(2,2,1) = CONJG(GCUSDL(1,1))
         GNFSF(1,2,1) = GNDSDL(1,1)
         GNFSF(2,2,1) = GNDSDL(2,1)
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      CALL OXXXXX(PV(0,1),AM1   ,IHEL(1),-1,SPINOR(1, 1))
      CALL IXXXXX(PV(0,2),AM2   ,IHEL(2),-1,SPINOR(1, 2))
      CALL OXXXXX(PV(0,3),AM3   ,IHEL(3),+1,SPINOR(1, 3))
      CALL IXXXXX(PV(0,3),AM3   ,IHEL(3),-1,SPINOR(1,13))
      CALL IXXXXX(PV(0,4),AM4   ,IHEL(4),-1,SPINOR(1, 4))
      CALL OXXXXX(PV(0,4),AM4   ,IHEL(4),+1,SPINOR(1,14))
C--
C  Prepare spinors for X+.
C--
      CALL JIOXXX(SPINOR(1,4),SPINOR(1,3),GWF,AMSW,GAMW,WVCT)
      CALL IOVXXC(SPINOR(1,2),SPINOR(1,1),WVCT,GWX(1,1),AMP(1))
C--
      CALL HIOXXX(SPINOR(1,2),SPINOR(1, 3),GNFSF(1,1,1),
     .            AMSSF(1,1),GAMSF(1,1),SCAL)
      CALL IOSXXX(SPINOR(1, 4),SPINOR(1,1),SCAL,GCFSF(1,1,1),AMP(2))
C--
      CALL HIOXXX(SPINOR(1,2),SPINOR(1,14),GNFSF(1,2,1),
     .            AMSSF(2,1),GAMSF(2,1),SCAL)
      CALL IOSXXX(SPINOR(1,13),SPINOR(1,1),SCAL,GCFSF(1,2,1),AMP(3))
C--
      AMP(0)  = AMP(1) - AMP(2) + AMP(3)
C>>>
CTBW       AMP(0) = AMP(1)
C>>>
C--
C  Calculate amplitude squared.
C--      
      FACT    = (3*QCD)**(IDP(3,3)-1)
     .          * VKM( IDP(1,3), IDP(1,4), IDP(3,3) )**2
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
      AMP2(3) = FACT*ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
