CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULCXX(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> XC+ + XC-.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  LSP from XC+.
C*                        :   =  4 for  fu from XC+
C*                        :   =  5 for  fdb from XC+.
C*                        :   =  6 for  LSP from XC-.
C*                        :   =  7 for  fub  from XC-
C*                        :   =  8 for  fd from XC-.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : photon.
C*        (2)     : (R*4) : Z.
C*        (3)     : (R*4) : t-channel snu-exchange.
C* (Relation)
C*    Invokes AMPCXX and subroutines in HELAS.LOAD.
C* (Update Record)
C*    92/06/09  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULCXX(IDP,IHEL,PV,AMP2)
 
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
      COMPLEX*8  SPINOR(6,20), EIN(6), EOT(6), XIN(6), XOT(6),
     .           WVCT(6), SCAL(3), AMP(0:3)
      COMPLEX*8  GWX(2,2), GSNX(2,2), GCFSF(2,2,2), GNFSF(2,2,2)
      REAL   *4  GAX(2), GZX(2), AMSSF(2,2), GAMSF(2,2)
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL     = 1
         AMSE      = AMF(1,2,1)
         AMSW      = AMW
         GAMW      = GMWTOT
         AMSZ      = AMZ
         GAMZ      = GMZTOT
C>>>
         QCD       =  1 + xALFS/xPI
C        QCD       =  1
C>>>
         GWX(1,1)  = CONJG(GNCW(1,1,1))
         GWX(2,1)  = CONJG(GNCW(2,1,1))
         GWX(1,2)  = GNCW(1,1,1)
         GWX(2,2)  = GNCW(2,1,1)
         GAX(1)    = GAL(1)
         GAX(2)    = GAL(2)
         GZX(1)    = GCCZ(1,1,1)
         GZX(2)    = GCCZ(2,1,1)
         GSNX(1,1) = GCESNL(1,1)
         GSNX(2,1) = GCESNL(2,1)
         GSNX(1,2) = CONJG(GCESNL(2,1))
         GSNX(2,2) = CONJG(GCESNL(1,1))
      ENDIF
C--
C  Set masses and couplings depending on final states.
C--
      AM4 = AMF(IDP(1,4),IDP(2,4),IDP(3,4))
      AM5 = AMF(IDP(1,5),IDP(2,5),IDP(3,5))
      AM7 = AMF(IDP(1,7),IDP(2,7),IDP(3,7))
      AM8 = AMF(IDP(1,8),IDP(2,8),IDP(3,8))
C--
      IF ( IDP(3,4).EQ.1 ) THEN
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
      IF ( IDP(3,7).EQ.1 ) THEN
         AMSSF(  1,2) = SFM(1)
         GAMSF(  1,2) = GMSF(1)
         GCFSF(1,1,2) = CONJG(GCESNL(2,1))
         GCFSF(2,1,2) = CONJG(GCESNL(1,1))
         GNFSF(1,1,2) = GNNSNL(1,1)
         GNFSF(2,1,2) = GNNSNL(2,1)
         AMSSF(  2,2) = SFM(2)
         GAMSF(  2,2) = GMSF(2)
         GCFSF(1,2,2) = GCNSEL(1,1)
         GCFSF(2,2,2) = GCNSEL(2,1)
         GNFSF(1,2,2) = CONJG(GNESEL(2,1))
         GNFSF(2,2,2) = CONJG(GNESEL(1,1))
      ELSE
         AMSSF(  1,2) = SFM(4)
         GAMSF(  1,2) = GMSF(4)
         GCFSF(1,1,2) = CONJG(GCDSUL(2,1))
         GCFSF(2,1,2) = CONJG(GCDSUL(1,1))
         GNFSF(1,1,2) = GNUSUL(1,1)
         GNFSF(2,1,2) = GNUSUL(2,1)
         AMSSF(  2,2) = SFM(6)
         GAMSF(  2,2) = GMSF(6)
         GCFSF(1,2,2) = GCUSDL(1,1)
         GCFSF(2,2,2) = GCUSDL(2,1)
         GNFSF(1,2,2) = CONJG(GNDSDL(2,1))
         GNFSF(2,2,2) = CONJG(GNDSDL(1,1))
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      CALL IXXXXX(PV(0,1),AMSE  ,IHEL(1),+1,EIN)
      CALL OXXXXX(PV(0,2),AMSE  ,IHEL(2),-1,EOT)
C--
      CALL IXXXXX(PV(0,3),SZM(1),IHEL(3),-1,SPINOR(1, 3))
      CALL OXXXXX(PV(0,4),AM4   ,IHEL(4),+1,SPINOR(1, 4))
      CALL IXXXXX(PV(0,4),AM4   ,IHEL(4),-1,SPINOR(1,14))
      CALL IXXXXX(PV(0,5),AM5   ,IHEL(5),-1,SPINOR(1, 5))
      CALL OXXXXX(PV(0,5),AM5   ,IHEL(5),+1,SPINOR(1,15))
C--
      CALL OXXXXX(PV(0,6),SZM(1),IHEL(6),+1,SPINOR(1, 6))
      CALL IXXXXX(PV(0,7),AM7   ,IHEL(7),-1,SPINOR(1, 7))
      CALL OXXXXX(PV(0,7),AM7   ,IHEL(7),+1,SPINOR(1,17))
      CALL OXXXXX(PV(0,8),AM8   ,IHEL(8),+1,SPINOR(1, 8))
      CALL IXXXXX(PV(0,8),AM8   ,IHEL(8),-1,SPINOR(1,18))
C--
C  Prepare spinors for X+.
C--
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,4),GWF,AMSW,GAMW,WVCT)
      CALL FVIXXC(SPINOR(1,3),WVCT,GWX(1,1),SWM(1),GMSW(1),XIN)
CNOXCTOW       XIN(1) = 0
CNOXCTOW       XIN(2) = 0
CNOXCTOW       XIN(3) = 0
CNOXCTOW       XIN(4) = 0
C--
      CALL HIOXXX(SPINOR(1,3),SPINOR(1, 4),GNFSF(1,1,1),
     .            AMSSF(1,1),GAMSF(1,1),SCAL)
      CALL FSIXXX(SPINOR(1, 5),SCAL,GCFSF(1,1,1),
     .            SWM(1),GMSW(1),SPINOR(1, 1))
CNOUSF       SPINOR(1,1) = 0
CNOUSF       SPINOR(2,1) = 0
CNOUSF       SPINOR(3,1) = 0
CNOUSF       SPINOR(4,1) = 0
C--
      CALL HIOXXX(SPINOR(1,3),SPINOR(1,15),GNFSF(1,2,1),
     .            AMSSF(2,1),GAMSF(2,1),SCAL)
      CALL FSIXXX(SPINOR(1,14),SCAL,GCFSF(1,2,1),
     .            SWM(1),GMSW(1),SPINOR(1,11))
CNODSF       SPINOR(1,11) = 0
CNODSF       SPINOR(2,11) = 0
CNODSF       SPINOR(3,11) = 0
CNODSF       SPINOR(4,11) = 0
C--
      XIN(1) = XIN(1) - SPINOR(1,1) + SPINOR(1,11)
      XIN(2) = XIN(2) - SPINOR(2,1) + SPINOR(2,11)
      XIN(3) = XIN(3) - SPINOR(3,1) + SPINOR(3,11)
      XIN(4) = XIN(4) - SPINOR(4,1) + SPINOR(4,11)
C--
C  Prepare spinors for X-.
C--
      CALL JIOXXX(SPINOR(1,7),SPINOR(1,8),GWF,AMSW,GAMW,WVCT)
      CALL FVOXXC(SPINOR(1,6),WVCT,GWX(1,2),SWM(1),GMSW(1),XOT)
CNOXCTOW       XOT(1) = 0
CNOXCTOW       XOT(2) = 0
CNOXCTOW       XOT(3) = 0
CNOXCTOW       XOT(4) = 0
C--
      CALL HIOXXX(SPINOR(1, 7),SPINOR(1,6),GNFSF(1,1,2),
     .            AMSSF(1,2),GAMSF(1,2),SCAL)
      CALL FSOXXX(SPINOR(1, 8),SCAL,GCFSF(1,1,2),
     .            SWM(1),GMSW(1),SPINOR(1, 2))
CNOUSF       SPINOR(1,2) = 0
CNOUSF       SPINOR(2,2) = 0
CNOUSF       SPINOR(3,2) = 0
CNOUSF       SPINOR(4,2) = 0
C--
      CALL HIOXXX(SPINOR(1,18),SPINOR(1,6),GNFSF(1,2,2),
     .            AMSSF(2,2),GAMSF(2,2),SCAL)
      CALL FSOXXX(SPINOR(1,17),SCAL,GCFSF(1,2,2),
     .            SWM(1),GMSW(1),SPINOR(1,12))
CNODSF       SPINOR(1,12) = 0
CNODSF       SPINOR(2,12) = 0
CNODSF       SPINOR(3,12) = 0
CNODSF       SPINOR(4,12) = 0
C--
      XOT(1) = XOT(1) - SPINOR(1,2) + SPINOR(1,12)
      XOT(2) = XOT(2) - SPINOR(2,2) + SPINOR(2,12)
      XOT(3) = XOT(3) - SPINOR(3,2) + SPINOR(3,12)
      XOT(4) = XOT(4) - SPINOR(4,2) + SPINOR(4,12)
C--
C  Calculate amplitude**2.
C--
      CALL AMPCXX (GAL,GAX,GZL,GZX,GSNX,SFM(1),GMSF(1),AMSZ,GAMZ,
     .             EIN,EOT,XIN,XOT,AMP)
      FACT    = (3*QCD)**(IDP(3,4)+IDP(3,7)-2)
     .          * VKM( IDP(1,4), IDP(1,5), IDP(3,4) )**2
     .          * VKM( IDP(1,7), IDP(1,8), IDP(3,7) )**2
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
      AMP2(3) = FACT*ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
