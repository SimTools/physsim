CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULNNZZ(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> nu nubar Z1 Z2.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  nu.
C*                        :   =  4 for  nubar.
C*                        :   =  5 for  Z1.
C*                        :   =  6 for  Z2.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : non-fusion.
C*        (2)     : (R*8) : fusion.
C*        (3)     : (R*8) : Higgs.
C* (Relation)
C*    Invokes AMPNNZZ and subroutines in HELAS.LOAD.
C* (Update Record)
C*    2005/03/06  K.Fujii        Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULNNZZ(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      PARAMETER  ( NP = 6 )
      INTEGER*4  IDP(3,NP), IHEL(NP)
      REAL   *8  PV(0:3,NP), AMP2(0:3)
C--
      COMPLEX*16 SPINOR(6,NP-2), VCT(6,2), AMP(0:3)
      INTEGER*4  NF(NP-2,2)
      DATA  NF  / +1, -1,   +1, -1,
     .            +1, -1,   -1, +1 /
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
         AMSH      = AMH
         GAMH      = GMHTOT
C--
         GWWH      = GW*AMSW
         GZZH      = GW*AMSW/xCOS2W
C--
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      DO 100 IP = 1, NP-2
         IF ( NF(IP,2).GT.0 ) THEN
            CALL IXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ELSE
            CALL OXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ENDIF
100   CONTINUE
C--
C  Prepare wave function for Z's.
C--
      CALL VXXXXX(PV(0,5),AMSZ,IHEL(5),+1,VCT(1,1))
      CALL VXXXXX(PV(0,6),AMSZ,IHEL(6),+1,VCT(1,2))
C--
C  Calculate amplitude**2.
C--
      CALL AMPNNZZ(GZL, GZN, GWF, GWWZ, GWWH, GZZH,
     .             AMSW,GAMW,AMSZ,GAMZ,AMSE, AMSH,GAMH,   
     .                         SPINOR(1,1),SPINOR(1,2),
     .                         SPINOR(1,3),SPINOR(1,4),
     .                                        VCT, AMP)
C--
      AMP2(0) = ABS(AMP(0))**2
      AMP2(1) = ABS(AMP(1))**2
      AMP2(2) = ABS(AMP(2))**2
      AMP2(3) = ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
