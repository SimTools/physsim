CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULNXX(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate full amplitude for e+ + e- ---> X0_i + X0_j
C*    followed by X0_i -> tau + stau
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  tau-/+ from X0_i
C*                        :   =  4 for stau+/- from X0_i
C*                        :   =  5 for  tau+/- from X0_j
C*                        :   =  6 for stau-/+ from X0_j
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : s-channel Z.
C*        (2)     : (R*8) : t/u-channel seL.
C*        (3)     : (R*8) : t/u-channel seR.
C* (Relation)
C*    Invokes AMPNXX and subroutines in HELAS.LOAD.
C* (Update Record)
C*    2010/10/07  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULNXX(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SSCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SMCUPL.inc"
#include "SSCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *8  PV(0:3,*), AMP2(0:3)
C--
      COMPLEX*16 SPINOR(6,20), EIN(6), EOT(6), 
     .           XIIN(6), XJIN(6), XIOT(6), XJOT(6),
     .           ZVCT(6), SCAL(3,2), AMP(0:3)
      COMPLEX*16 GZX(2,6), GSEX(2,2,2,2), GSFX(2,2,2,2)
      REAL   *8  AMSSE(2), GAMSE(2), AMSSF(2), GAMSF(2)
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL     = 1
C>>>
         QCD       =  1 + xALFS/xPI
C        QCD       =  1
C>>>
C--
C  Masses and widths.
C--
         AMSE      = AMF(1,2,1)
         AMSZ      = AMZ
         GAMZ      = GMZTOT
C
C-- (1,2)=(seL,seR)
C
         AMSSE(1)  = SFM (2)
         GAMSE(1)  = GMSF(2)
         AMSSE(2)  = SFM (3)
         GAMSE(2)  = GMSF(3)
C--
C  SUSY couplings.
C     GZX (*,1)        : X0_J-X0_I-Z coupling.
C         (*,2)        : X0_I-X0_J-Z coupling.
C     GSEX(*,IX,LR,IO) : IO=(1,2)=(se_LR-e-X0_[IX], X0_[IX]-e-se_LR)
C                      : LR=(1,2)=(left,right)
C--
#if 0
         IX        = 2
         JX        = 2
#else
         IX        = 1
         JX        = 1
#endif
         GZX(1,1)       = -GNNZ(2,IX,JX)
         GZX(2,1)       = -GNNZ(1,IX,JX)
         GZX(1,2)       =  GNNZ(1,IX,JX)
         GZX(2,2)       =  GNNZ(2,IX,JX)
         GZX(1,3)       = -GNNZ(2,1,IX)
         GZX(2,3)       = -GNNZ(1,1,IX)
         GZX(1,4)       =  GNNZ(1,1,IX)
         GZX(2,4)       =  GNNZ(2,1,IX)
         GZX(1,5)       = -GNNZ(2,1,JX)
         GZX(2,5)       = -GNNZ(1,1,JX)
         GZX(1,6)       =  GNNZ(1,1,JX)
         GZX(2,6)       =  GNNZ(2,1,JX)
C--
         LR = 1
         GSEX(1,1,LR,2) = GNESEL(1,IX)
         GSEX(2,1,LR,2) = GNESEL(2,IX)
         GSEX(1,1,LR,1) = CONJG(GNESEL(2,IX))
         GSEX(2,1,LR,1) = CONJG(GNESEL(1,IX))
C--
         LR = 2
         GSEX(1,1,LR,2) = GNESER(1,IX)
         GSEX(2,1,LR,2) = GNESER(2,IX)
         GSEX(1,1,LR,1) = CONJG(GNESER(2,IX))
         GSEX(2,1,LR,1) = CONJG(GNESER(1,IX))
C--         
         LR = 1
         GSEX(1,2,LR,2) = GNESEL(1,JX)
         GSEX(2,2,LR,2) = GNESEL(2,JX)
         GSEX(1,2,LR,1) = CONJG(GNESEL(2,JX))
         GSEX(2,2,LR,1) = CONJG(GNESEL(1,JX))
C--
         LR = 2
         GSEX(1,2,LR,2) = GNESER(1,JX)
         GSEX(2,2,LR,2) = GNESER(2,JX)
         GSEX(1,2,LR,1) = CONJG(GNESER(2,JX))
         GSEX(2,2,LR,1) = CONJG(GNESER(1,JX))
C--
         PRINT *, '*************************************'
         PRINT *, ' AMX(1)   = ', SZM(IX)
         PRINT *, '    (2)   = ', SZM(JX)
         PRINT *, ' GMX(1)   = ', GMSZ(IX)
         PRINT *, '    (2)   = ', GMSZ(JX)
         PRINT *, ' AMSZ     = ', AMSZ
         PRINT *, ' GAMZ     = ', GAMZ
         PRINT *, ' AMSSE(1) = ', AMSSE(1)
         PRINT *, '      (2) = ', AMSSE(2) 
         PRINT *, ' GAMSE(1) = ', GAMSE(1)
         PRINT *, '      (2) = ', GAMSE(2) 
         PRINT *, ' '
         PRINT *, ' GZL(1)   = ', GZL(1)
         PRINT *, '    (2)   = ', GZL(2)
         PRINT *, ' GZX(1,1) = ', GZX(1,1)
         PRINT *, '    (2,1) = ', GZX(2,1)
         PRINT *, ' GZX(1,2) = ', GZX(1,2)
         PRINT *, '    (2,2) = ', GZX(2,2)
         PRINT *, ' GZX(1,3) = ', GZX(1,3)
         PRINT *, '    (2,3) = ', GZX(2,3)
         PRINT *, ' GZX(1,4) = ', GZX(1,4)
         PRINT *, '    (2,4) = ', GZX(2,4)
         PRINT *, ' GZX(1,5) = ', GZX(1,3)
         PRINT *, '    (2,5) = ', GZX(2,3)
         PRINT *, ' GZX(1,6) = ', GZX(1,4)
         PRINT *, '    (2,6) = ', GZX(2,4)
         PRINT *, ' '
         PRINT *, ' GSEX(1,1,1,1) = ', GSEX(1,1,1,1)
         PRINT *, '     (2,1,1,1) = ', GSEX(2,1,1,1)
         PRINT *, '     (1,2,1,1) = ', GSEX(1,2,1,1)
         PRINT *, '     (2,2,1,1) = ', GSEX(2,2,1,1)
         PRINT *, ' GSEX(1,1,2,1) = ', GSEX(1,1,2,1)
         PRINT *, '     (2,1,2,1) = ', GSEX(2,1,2,1)
         PRINT *, '     (1,2,2,1) = ', GSEX(1,2,2,1)
         PRINT *, '     (2,2,2,1) = ', GSEX(2,2,2,1)
         PRINT *, ' GSEX(1,1,1,2) = ', GSEX(1,1,1,2)
         PRINT *, '     (2,1,1,2) = ', GSEX(2,1,1,2)
         PRINT *, '     (1,2,1,2) = ', GSEX(1,2,1,2)
         PRINT *, '     (2,2,1,2) = ', GSEX(2,2,1,2)
         PRINT *, ' GSEX(1,1,2,2) = ', GSEX(1,1,2,2)
         PRINT *, '     (2,1,2,2) = ', GSEX(2,1,2,2)
         PRINT *, '     (1,2,2,2) = ', GSEX(1,2,2,2)
         PRINT *, '     (2,2,2,2) = ', GSEX(2,2,2,2)
         PRINT *, '*************************************'
      ENDIF
C--
C  Set masses and couplings depending on final states.
C  This version ignores E_L, U_L, D_L, U_2, and D_2,
C  assuming that they are heavier than X0_j.
C--
C  SUSY couplings.
C     GSFX(*,NP,I12,IO) : IO=(1,2)=(sf-t-X0_[IX], X0_[IX]-t-sf)
C                       : IJ=(1,2)=(I,J)
C                       : NP=(1,2)=(particle,anti-particle)
C--
C  SFM(i) : i=(1,2,3,4,5,6,7) = (snl,sel,ser,sul,sur,sdl,sdr)
C--
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
               AMSSF(1) = SFM3M(2)
               GAMSF(1) = GMSF3(2)
C--
               GSFX(1,1,1,2) = GNESE1(1,IX)
               GSFX(2,1,1,2) = GNESE1(2,IX)
               GSFX(1,1,1,1) = CONJG(GNESE1(2,IX))
               GSFX(2,1,1,1) = CONJG(GNESE1(1,IX))
C--
               GSFX(1,2,1,2) = CONJG(GNESE1(2,IX))
               GSFX(2,2,1,2) = CONJG(GNESE1(1,IX))
               GSFX(1,2,1,1) = GNESE1(1,IX)
               GSFX(2,2,1,1) = GNESE1(2,IX)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C--
      AM3 = AMF(IDP(1,3),IDP(2,3),IDP(3,3))
      AM4 = AMSSF(1)
C--
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
               AMSSF(2) = SFM3M(2)
               GAMSF(2) = GMSF3(2)
C--
               GSFX(1,1,2,2) = GNESE1(1,JX)
               GSFX(2,1,2,2) = GNESE1(2,JX)
               GSFX(1,1,2,1) = CONJG(GNESE1(2,JX))
               GSFX(2,1,2,1) = CONJG(GNESE1(1,JX))
C--
               GSFX(1,2,2,2) = CONJG(GNESE1(2,JX))
               GSFX(2,2,2,2) = CONJG(GNESE1(1,JX))
               GSFX(1,2,2,1) = GNESE1(1,JX)
               GSFX(2,2,2,1) = GNESE1(2,JX)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C--
      AM5 = AMF(IDP(1,5),IDP(2,5),IDP(3,5))
      AM6 = AMSSF(2)
C--
C  Prepare spinors for external lines.
C--
      CALL IXXXXX(PV(0,1),AMSE  ,IHEL(1),+1,EIN)
      CALL OXXXXX(PV(0,2),AMSE  ,IHEL(2),-1,EOT)
C-- 
      IF (IDP(2,4).GT.0) THEN
         CALL OXXXXX(PV(0,3),AM3   ,IHEL(3),+1,SPINOR(1, 3))
         CALL IXXXXX(PV(0,3),AM3   ,IHEL(3),-1,SPINOR(1,13))
      ELSE
         CALL IXXXXX(PV(0,3),AM3   ,IHEL(3),-1,SPINOR(1, 3))
         CALL OXXXXX(PV(0,3),AM3   ,IHEL(3),+1,SPINOR(1,13))
      ENDIF
      CALL SXXXXX(PV(0,4), +1   ,           SCAL  (1, 1))

      IF (IDP(2,6).GT.0) THEN
         CALL IXXXXX(PV(0,5),AM5   ,IHEL(5),-1,SPINOR(1, 5))
         CALL OXXXXX(PV(0,5),AM5   ,IHEL(5),+1,SPINOR(1,15))
      ELSE
         CALL OXXXXX(PV(0,5),AM5   ,IHEL(5),+1,SPINOR(1, 5))
         CALL IXXXXX(PV(0,5),AM5   ,IHEL(5),-1,SPINOR(1,15))
      ENDIF
      CALL SXXXXX(PV(0,6), +1   ,           SCAL  (1, 2))
C--
C  Prepare spinors for X0_IX.
C--
C                                             *,np,IX/JX,io
      IF (IDP(2,4).GT.0) THEN
         CALL FSOXXX(SPINOR(1, 3),SCAL(1,1),GSFX(1,1,1,1),
     .               SZM(IX),GMSZ(IX), XIOT)
         CALL FSIXXX(SPINOR(1,13),SCAL(1,1),GSFX(1,2,1,2),
     .               SZM(IX),GMSZ(IX), XIIN)
      ELSE
         CALL FSIXXX(SPINOR(1, 3),SCAL(1,1),GSFX(1,1,1,2),
     .               SZM(IX),GMSZ(IX), XIIN)
         CALL FSOXXX(SPINOR(1,13),SCAL(1,1),GSFX(1,2,1,1),
     .               SZM(IX),GMSZ(IX), XIOT)
      ENDIF
C--
C  Prepare spinors for X0_JX.
C--
      IF (IDP(2,6).GT.0) THEN
         CALL FSIXXX(SPINOR(1, 5),SCAL(1,2),GSFX(1,1,2,2),
     .               SZM(JX),GMSZ(JX), XJIN)
         CALL FSOXXX(SPINOR(1,15),SCAL(1,2),GSFX(1,2,2,1),
     .               SZM(JX),GMSZ(JX), XJOT)
      ELSE
         CALL FSOXXX(SPINOR(1, 5),SCAL(1,2),GSFX(1,1,2,1),
     .               SZM(JX),GMSZ(JX), XJOT)
         CALL FSIXXX(SPINOR(1,15),SCAL(1,2),GSFX(1,2,2,2),
     .               SZM(JX),GMSZ(JX), XJIN)
      ENDIF
C--
C  Calculate amplitude**2.
C--
CCDBG        PRINT *, ' EIN = ', EIN
CCDBG        PRINT *, ' EOT = ', EOT
CCDBG        PRINT *, ' XIIN = ', XIIN
CCDBG        PRINT *, ' XIOT = ', XIOT
CCDBG        PRINT *, ' XJIN = ', XJIN
CCDBG        PRINT *, ' XJOT = ', XJOT
      CALL AMPNXX (GZL,GZX(1,2),GSEX,AMSSE,GAMSE,AMSZ,GAMZ,
     .             EIN,EOT,XIIN,XJIN,XIOT,XJOT, AMP)
CCDBG        PRINT *, ' AMP = ', AMP
      FACT    = (3*QCD)**(IDP(3,3)+IDP(3,5)-2)
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
      AMP2(3) = FACT*ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
