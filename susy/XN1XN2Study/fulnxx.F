CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULNXX(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate full amplitude for e+ + e- ---> X0_1 + X0_2
C*    followed by X0_2 -> f + fbar + X0_1 
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  X0_1 = LSP.
C*                        :   =  4 for  tau- from X0_2
C*                        :   =  5 for  tau+ from X0_2
C*                        :   =  6 for  LSP  from X0=2
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : s-channel Z.
C*        (2)     : (R*8) : t/u-channel seL.
C*        (3)     : (R*8) : t/u-channel seR.
C* (Relation)
C*    Invokes AMPNXX and subroutines in HELAS.LOAD.
C* (Update Record)
C*    03/08/28  K.Fujii   Original version.
C*    03/09/10  K.Fujii	  Modified to allow
C*                        Z2->snL+nu, seR+e, suR+u, sdR+d (1st and 2nd)
C*                          ->snL+nu, se1+e, su1+u, sd1+d (3rd gen)
C*                        !!! CAUTION !!!
C*                          seL, se2, suL, sdL decays are ignored.
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
      COMPLEX*16 GZX(2,2), GSEX(2,2,2,2), GSFX(2,2,2)
      REAL   *8  GZF(2)
      REAL   *8  AMSSE(2), GAMSE(2), AMSSF, GAMSF
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
C--
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
         IX        = 1
         JX        = 2
         GZX(1,1)   = -GNNZ(2,IX,JX)
         GZX(2,1)   = -GNNZ(1,IX,JX)
         GZX(1,2)   = GNNZ(1,IX,JX)
         GZX(2,2)   = GNNZ(2,IX,JX)
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
      ENDIF
C--
C  Set masses and couplings depending on final states.
C  This version ignores E_L, U_L, D_L, U_2, and D_2,
C  assuming that they are heavier than X0_j.
C--
C  SUSY couplings.
C     GSFX(*,IX,IO) : IO=(1,2)=(sf-t-X0_[IX], X0_[IX]-t-sf)
C                   : IX=(1,2)=(I,J)
C--
      IF (IDP(3,4).EQ.1) THEN
         IF (IDP(1,4).NE.3) THEN
            IF (IDP(2,4).EQ.1) THEN
               AMSSF  = SFM (1)
               GAMSF  = GMSF(1)
C--
               GZF(1) = GZN(1)
               GZF(2) = GZN(2)
C--
               GSFX(1,1,2) = GNNSNL(1,IX)
               GSFX(2,1,2) = GNNSNL(2,IX)
               GSFX(1,1,1) = CONJG(GNNSNL(2,IX))
               GSFX(2,1,1) = CONJG(GNNSNL(1,IX))
C--         
               GSFX(1,2,2) = GNNSNL(1,JX)
               GSFX(2,2,2) = GNNSNL(2,JX)
               GSFX(1,2,1) = CONJG(GNNSNL(2,JX))
               GSFX(2,2,1) = CONJG(GNNSNL(1,JX))
            ELSE
               AMSSF  = SFM (3)
               GAMSF  = GMSF(3)
C--
               GZF(1) = GZL(1)
               GZF(2) = GZL(2)
C--
               GSFX(1,1,2) = GNESER(1,IX)
               GSFX(2,1,2) = GNESER(2,IX)
               GSFX(1,1,1) = CONJG(GNESER(2,IX))
               GSFX(2,1,1) = CONJG(GNESER(1,IX))
C--
               GSFX(1,2,2) = GNESER(1,JX)
               GSFX(2,2,2) = GNESER(2,JX)
               GSFX(1,2,1) = CONJG(GNESER(2,JX))
               GSFX(2,2,1) = CONJG(GNESER(1,JX))
            ENDIF
         ELSE
            IF (IDP(2,4).EQ.1) THEN
               AMSSF  = SFM3M(1)
               GAMSF  = GMSF3(1)
C--
               GZF(1) = GZN(1)
               GZF(2) = GZN(2)
C--
               GSFX(1,1,2) = GNSNL3(1,IX)
               GSFX(2,1,2) = GNSNL3(2,IX)
               GSFX(1,1,1) = CONJG(GNSNL3(2,IX))
               GSFX(2,1,1) = CONJG(GNSNL3(1,IX))
C--         
               GSFX(1,2,2) = GNSNL3(1,JX)
               GSFX(2,2,2) = GNSNL3(2,JX)
               GSFX(1,2,1) = CONJG(GNSNL3(2,JX))
               GSFX(2,2,1) = CONJG(GNSNL3(1,JX))
            ELSE
               AMSSF  = SFM3M(2)
               GAMSF  = GMSF3(2)
C--
               GZF(1) = GZL(1)
               GZF(2) = GZL(2)
C--
               GSFX(1,1,2) = GNESE1(1,IX)
               GSFX(2,1,2) = GNESE1(2,IX)
               GSFX(1,1,1) = CONJG(GNESE1(2,IX))
               GSFX(2,1,1) = CONJG(GNESE1(1,IX))
C--
               GSFX(1,2,2) = GNESE1(1,JX)
               GSFX(2,2,2) = GNESE1(2,JX)
               GSFX(1,2,1) = CONJG(GNESE1(2,JX))
               GSFX(2,2,1) = CONJG(GNESE1(1,JX))
            ENDIF
         ENDIF
      ELSE
         IF (IDP(1,4).NE.3) THEN
            IF (IDP(2,4).EQ.1) THEN
               AMSSF  = SFM (5)
               GAMSF  = GMSF(5)
C--
               GZF(1) = GZU(1)
               GZF(2) = GZU(2)
C--
               GSFX(1,1,2) = GNUSUR(1,IX)
               GSFX(2,1,2) = GNUSUR(2,IX)
               GSFX(1,1,1) = CONJG(GNUSUR(2,IX))
               GSFX(2,1,1) = CONJG(GNUSUR(1,IX))
C--         
               GSFX(1,2,2) = GNUSUR(1,JX)
               GSFX(2,2,2) = GNUSUR(2,JX)
               GSFX(1,2,1) = CONJG(GNUSUR(2,JX))
               GSFX(2,2,1) = CONJG(GNUSUR(1,JX))
            ELSE
               AMSSF  = SFM (7)
               GAMSF  = GMSF(7)
C--
               GZF(1) = GZD(1)
               GZF(2) = GZD(2)
C--
               GSFX(1,1,2) = GNDSDR(1,IX)
               GSFX(2,1,2) = GNDSDR(2,IX)
               GSFX(1,1,1) = CONJG(GNDSDR(2,IX))
               GSFX(2,1,1) = CONJG(GNDSDR(1,IX))
C--
               GSFX(1,2,2) = GNDSDR(1,JX)
               GSFX(2,2,2) = GNDSDR(2,JX)
               GSFX(1,2,1) = CONJG(GNDSDR(2,JX))
               GSFX(2,2,1) = CONJG(GNDSDR(1,JX))
            ENDIF
         ELSE
            IF (IDP(2,4).EQ.1) THEN
               AMSSF  = SFM3M(4)
               GAMSF  = GMSF3(4)
C--
               GZF(1) = GZU(1)
               GZF(2) = GZU(2)
C--
               GSFX(1,1,2) = GNUSU1(1,IX)
               GSFX(2,1,2) = GNUSU1(2,IX)
               GSFX(1,1,1) = CONJG(GNUSU1(2,IX))
               GSFX(2,1,1) = CONJG(GNUSU1(1,IX))
C--         
               GSFX(1,2,2) = GNUSU1(1,JX)
               GSFX(2,2,2) = GNUSU1(2,JX)
               GSFX(1,2,1) = CONJG(GNUSU1(2,JX))
               GSFX(2,2,1) = CONJG(GNUSU1(1,JX))
            ELSE
               AMSSF  = SFM3M(6)
               GAMSF  = GMSF3(6)
C--
               GZF(1) = GZD(1)
               GZF(2) = GZD(2)
C--
               GSFX(1,1,2) = GNDSD1(1,IX)
               GSFX(2,1,2) = GNDSD1(2,IX)
               GSFX(1,1,1) = CONJG(GNDSD1(2,IX))
               GSFX(2,1,1) = CONJG(GNDSD1(1,IX))
C--
               GSFX(1,2,2) = GNDSD1(1,JX)
               GSFX(2,2,2) = GNDSD1(2,JX)
               GSFX(1,2,1) = CONJG(GNDSD1(2,JX))
               GSFX(2,2,1) = CONJG(GNDSD1(1,JX))
            ENDIF
         ENDIF
      ENDIF
      AM4 = AMF(IDP(1,4),IDP(2,4),IDP(3,4))
      AM5 = AMF(IDP(1,5),IDP(2,5),IDP(3,5))
C--
C  Prepare spinors for external lines.
C--
      CALL IXXXXX(PV(0,1),AMSE  ,IHEL(1),+1,EIN)
      CALL OXXXXX(PV(0,2),AMSE  ,IHEL(2),-1,EOT)
C--
      CALL IXXXXX(PV(0,3),SZM(1),IHEL(3),-1,XIIN)
      CALL OXXXXX(PV(0,3),SZM(1),IHEL(3),+1,XIOT)
C-- 
      CALL OXXXXX(PV(0,4),AM4   ,IHEL(4),+1,SPINOR(1, 4))
      CALL IXXXXX(PV(0,4),AM4   ,IHEL(4),-1,SPINOR(1,14))
      CALL IXXXXX(PV(0,5),AM5   ,IHEL(5),-1,SPINOR(1, 5))
      CALL OXXXXX(PV(0,5),AM5   ,IHEL(5),+1,SPINOR(1,15))
      CALL OXXXXX(PV(0,6),SZM(1),IHEL(6),+1,SPINOR(1, 6))
      CALL IXXXXX(PV(0,6),SZM(1),IHEL(6),-1,SPINOR(1,16))
C--
C  Prepare wave functions for Z, stau1-, and stau1+.
C--
      CALL JIOXXX(SPINOR(1, 5),SPINOR(1, 4),GZF,AMSZ,GAMZ,ZVCT)
C-- sf
      CALL HIOXXX(SPINOR(1,16),SPINOR(1, 4),GSFX(1,1,1),
     .            AMSSF,GAMSF, SCAL(1,1))
C-- sfbar
      CALL HIOXXX(SPINOR(1, 5),SPINOR(1, 6),GSFX(1,1,2),
     .            AMSSF,GAMSF, SCAL(1,2))
C--
C  Prepare spinors for incoming X0_2.
C--
      CALL FVIXXC(SPINOR(1,16),ZVCT,GZX(1,1),
     .            SZM(JX),GMSZ(JX), XJIN)
C--
      CALL FSIXXX(SPINOR(1, 5),SCAL(1,1),GSFX(1,JX,2),
     .            SZM(JX),GMSZ(JX), SPINOR(1,11))
C--
      CALL FSIXXX(SPINOR(1,14),SCAL(1,2),GSFX(1,JX,1),
     .            SZM(JX),GMSZ(JX), SPINOR(1,12))
C--
      XJIN(1) = XJIN(1) - SPINOR(1,11) - SPINOR(1,12)
      XJIN(2) = XJIN(2) - SPINOR(2,11) - SPINOR(2,12)
      XJIN(3) = XJIN(3) - SPINOR(3,11) - SPINOR(3,12)
      XJIN(4) = XJIN(4) - SPINOR(4,11) - SPINOR(4,12)
C--
C  Prepare spinors for outgoing X0_2.
C--
      CALL FVOXXC(SPINOR(1, 6),ZVCT,GZX(1,2),
     .            SZM(JX),GMSZ(JX), XJOT)
C--
      CALL FSOXXX(SPINOR(1,15),SCAL(1,1),GSFX(1,JX,2),
     .            SZM(JX),GMSZ(JX), SPINOR(1, 1))
C--
      CALL FSOXXX(SPINOR(1, 4),SCAL(1,2),GSFX(1,JX,1),
     .            SZM(JX),GMSZ(JX), SPINOR(1, 2))
C--
      XJOT(1) = XJOT(1) - SPINOR(1, 1) - SPINOR(1, 2)
      XJOT(2) = XJOT(2) - SPINOR(2, 1) - SPINOR(2, 2)
      XJOT(3) = XJOT(3) - SPINOR(3, 1) - SPINOR(3, 2)
      XJOT(4) = XJOT(4) - SPINOR(4, 1) - SPINOR(4, 2)
C--
C  Calculate amplitude**2.
C--
CCDBG        PRINT *, ' EIN = ', EIN
CCDBG        PRINT *, ' EOT = ', EOT
CCDBG        PRINT *, ' XIIN = ', XIIN
CCDBG        PRINT *, ' XIOT = ', XIOT
CCDBG        PRINT *, ' XJIN = ', XJIN
CCDBG        PRINT *, ' XJOT = ', XJOT
      CALL AMPNXX (GZL,GZX,GSEX,AMSSE,GAMSE,AMSZ,GAMZ,
     .             EIN,EOT,XIIN,XJIN,XIOT,XJOT, AMP)
CCDBG        PRINT *, ' AMP = ', AMP
      FACT    = (3*QCD)**(IDP(3,4)-1)
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
      AMP2(3) = FACT*ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
