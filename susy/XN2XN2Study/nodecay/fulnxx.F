CC**********************************************************************
C*
C*=========================----======----===
C* Subroutine FULNXX(IX1,IX2,IHEL,PV,AMP2)
C*=========================----======----===
C*
C* (Purpose)
C*    Calculate full amplitude for e+ + e- ---> X0_i + X0_j.
C* (Inputs)
C*       IX1      : (I*4) : i
C*       IX2      : (I*4) : j
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  X0_i
C*                        :   =  4 for  X0_j
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : s-channel Z.
C*        (2)     : (R*8) : t/u-channel seL.
C*        (3)     : (R*8) : t/u-channel seR.
C* (Relation)
C*    Invokes AMPNXX and subroutines in HELAS.LOAD.
C* (Update Record)
C*    03/09/03  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULNXX(IX1,IX2,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SSCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SMCUPL.inc"
#include "SSCUPL.inc"
      INTEGER*4  IX1, IX2, IHEL(*)
      REAL   *8  PV(0:3,*), AMP2(0:3)
C--
      COMPLEX*16 EIN(6), EOT(6), XIIN(6), XJIN(6), XIOT(6), XJOT(6),
     .           AMP(0:3)
      COMPLEX*16 GZX(2,4), GSEX(2,2,2,2)
      REAL   *8  AMSSE(2), GAMSE(2)
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL     = 1
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
         IX        = 2
         JX        = 2
         GZX(1,1)       = -GNNZ(2,IX,JX)
         GZX(2,1)       = -GNNZ(1,IX,JX)
         GZX(1,2)       =  GNNZ(1,IX,JX)
         GZX(2,2)       =  GNNZ(2,IX,JX)
         GZX(1,3)       = -GNNZ(2,1,JX)
         GZX(2,3)       = -GNNZ(1,1,JX)
         GZX(1,4)       =  GNNZ(1,1,JX)
         GZX(2,4)       =  GNNZ(2,1,JX)
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
C  Prepare spinors for external lines.
C--
      CALL IXXXXX(PV(0,1),AMSE  ,IHEL(1),+1,EIN)
      CALL OXXXXX(PV(0,2),AMSE  ,IHEL(2),-1,EOT)
C-- 
      CALL IXXXXX(PV(0,3),SZM(IX1),IHEL(3),-1,XIIN)
      CALL OXXXXX(PV(0,3),SZM(IX1),IHEL(3),+1,XIOT)
C--
      CALL IXXXXX(PV(0,4),SZM(IX2),IHEL(4),-1,XJIN)
      CALL OXXXXX(PV(0,4),SZM(IX2),IHEL(4),+1,XJOT)
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
      AMP2(0) = ABS(AMP(0))**2
      AMP2(1) = ABS(AMP(1))**2
      AMP2(2) = ABS(AMP(2))**2
      AMP2(3) = ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
