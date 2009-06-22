CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULZZZ(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> Z1 + Z2 + Z with
C*    decays of Z1, Z2, and Z into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  f from Z1.
C*                        :   =  4 for  f-bar from Z1.
C*                        :   =  5 for  f-bar from Z2.
C*                        :   =  6 for  f from Z2.
C*                        :   =  7 for  f-bar from Z.
C*                        :   =  8 for  f from Z.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : amp**2 of non-Higgs diagrams.
C*        (2)     : (R*4) : amp**2 of Higgs diagram.
C* (Relation)
C*    Invokes AMPZZZ and subroutines in HELAS.LOAD.
C* (Update Record)
C*    09/06/18  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULZZZ(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,8), IHEL(8)
      REAL   *4  PV(0:3,8), AMP2(0:2)
C--
      REAL   *4  GZF(2,2,2)
      COMPLEX*8  SPINOR(6,8), VCT(6,3), AMP(0:2)
      INTEGER*4  NF(8,2)
      DATA  NF  / +1, -1,   +1, -1,   -1, +1,   -1, +1,
     .            +1, -1,   -1, +1,   +1, -1,   +1, -1 /
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         CS2W   = xCOS2W
         AMSE   = AMF(1,2,1)
         AMSW   = AMW
         GAMW   = GMWTOT
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         AMSH   = AMH
         GAMH   = GMHTOT
C--
         GZH    = GW*AMSW/CS2W
         GZF(1,1,1) = GZN(1)
         GZF(2,1,1) = GZN(2)
         GZF(1,2,1) = GZL(1)
         GZF(2,2,1) = GZL(2)
         GZF(1,1,2) = GZU(1)
         GZF(2,1,2) = GZU(2)
         GZF(1,2,2) = GZD(1)
         GZF(2,2,2) = GZD(2)
         COLOR      = 3*(1+xALFS/xPI)
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      DO 100 IP = 1, 8
         IF ( NF(IP,2).GT.0 ) THEN
            CALL IXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ELSE
            CALL OXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ENDIF
100   CONTINUE
C--
C  Prepare wave function for Z1 and Z2.
C--
      CALL JIOXXX(SPINOR(1,4),SPINOR(1,3),
     .            GZF(1,IDP(2,4),IDP(3,4)),AMSZ,GAMZ,VCT(1,1))
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,6),
     .            GZF(1,IDP(2,5),IDP(3,5)),AMSZ,GAMZ,VCT(1,2))
C--
C  Prepare wave function for Z0.
C--
      CALL JIOXXX(SPINOR(1,7),SPINOR(1,8),
     .            GZF(1,IDP(2,7),IDP(3,7)),AMSZ,GAMZ,VCT(1,3))
C--
C  Calculate amplitude**2.
C--
      CALL AMPZZZ(GAL,GZN,GZL,GWF,GZH,GWH,GW,GZZZ,
     .            AMSE,AMSW,GAMW,AMSZ,GAMZ,AMSH,GAMH,
     .            VCT, SPINOR(1,1),SPINOR(1,2), AMP)
      FACT    = COLOR**(IDP(3,3)+IDP(3,5)+IDP(3,7)-3)
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END
