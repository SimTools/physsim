CC**********************************************************************
C*
C*==============================----===
C* Subroutine FULTBW(IDP,IHEL,PV,AMP2)
C*==============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + b-bar * W- with
C*    decays of t and t-bar into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  bbar.
C*                        :   =  4 for  t.
C*                        :   =  5 for  W-.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : tt
C*        (2)     : (R*8) : bb 
C*        (3)     : (R*8) : ww
C* (Relation)
C*    Invokes AMPTBW and subroutines in HELAS.LOAD.
C* (Update Record)
C*   2010/11/29  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULTBW(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *8  PV(0:3,*), AMP2(0:3)
C--
      REAL   *8  GAL(2), GAU(2), GAD(2), GWF(2), GZN(2), GZL(2), GZU(2),
     .           GZD(2), G1(2), GW, GWWZ, GWWA
      COMPLEX*16 SPINOR(6,4), WM(6), AMP(0:3)
      INTEGER*4  NF(4,2)
      DATA  NF  / +1, -1,   -1, +1,
     .            +1, -1,   +1, -1/
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         ALF    = xALF
         SN2W   = xSIN2W
         CS2W   = xCOS2W
         AMSE   = AMF(1,2,1)
         AMST   = AMF(3,1,2)
         GAMT   = GMF(3,1,2)
         AMSB   = AMF(3,2,2)
         GAMB   = GMF(3,2,2)
         AMSH   = AMH
         GAMH   = GMHTOT
         AMSW   = AMW
         GAMW   = GMWTOT
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         CALL UHCUPL(ALF,SN2W, GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1,
     .               GW,GWWZ,GWWA)
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      DO 100 IP = 1, 4
         IF ( NF(IP,2).GT.0 ) THEN
            CALL IXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ELSE
            CALL OXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ENDIF
100   CONTINUE
C-- W-
      CALL VXXXXX(PV(0,5),AMSW,IHEL(5),+1,WM)
C--
C  Calculate amplitude**2.
C  AMPTBW outputs
C     AMP(0) = AMP(1) + AMP(2) + AMP(3)
C        (1) = tt part
C        (2) = bb part.
C        (3) = ww part.
C--
      CALL AMPTBW(GAL,GAU,GAD,GZL,GZU,GZD,GWF,GW,AMSW,GAMW,AMSZ,GAMZ,
     .            AMST,GAMT,AMSB,GAMB,
     .            SPINOR(1, 1), SPINOR(1, 2),
     .            SPINOR(1, 3), SPINOR(1, 4), WM, AMP)
      FACT    = 3
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
      AMP2(3) = FACT*ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END

