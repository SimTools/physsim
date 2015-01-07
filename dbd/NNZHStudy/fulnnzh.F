CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULNNZH(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> nu + nubar + Z + H.
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
C*                        :   =  5 for  fdb from Z.
C*                        :   =  6 for  fu from Z.
C*                        :   =  7 for  H.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : nnZ* followed by Z* -> ZH
C*        (2)     : (R*8) : sum of diagrams, where Z from e-,e+,nu,nubar
C*        (3)     : (R*8) : sum of diagrams, where Z from W near e-, W near e+
C* (Relation)
C*    Invokes AMPNNZH and subroutines in HELAS.LOAD.
C* (Update Record)
C*    2010/06/19  J.Tian & K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULNNZH(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *8  PV(0:3,*), AMP2(0:3)
C--
      COMPLEX*16 SPINOR(6,6), ZVCT(6), SCH(3), AMP(0:7)
      REAL   *8  GZF(2,2,2)
      INTEGER*4  NF(6,2)
      DATA  NF  / +1, -1,   +1, -1,   -1, +1,
     .            +1, -1,   -1, +1,   +1, -1 /
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
C--
         GWWH      = GW*AMSW
         GZZH      = GW*AMSW/xCOS2W
C--
C        GZF(L/R,up/down,lepton/quark)
         GZF(1,1,1) = GZN(1)
         GZF(2,1,1) = GZN(2)
         GZF(1,2,1) = GZL(1)
         GZF(2,2,1) = GZL(2)
         GZF(1,1,2) = GZU(1)
         GZF(2,1,2) = GZU(2)
         GZF(1,2,2) = GZD(1)
         GZF(2,2,2) = GZD(2)
C--
         COLOR     = 3
         QCD       = 1 + xALFS/xPI
C--
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      DO 100 IP = 1, 6
         IF ( NF(IP,2).GT.0 ) THEN
            CALL IXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ELSE
            CALL OXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ENDIF
100   CONTINUE
C--
C  Prepare wave function for Z.
C--
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,6), 
     .            GZF(1,IDP(2,5),IDP(3,5)),AMSZ,GAMZ, ZVCT(1))
C--
C  Prepare wave function for H.
C--
      CALL SXXXXX(PV(0,7),+1,SCH)
C--
C  Calculate amplitude**2.
C--
      CALL AMPNNZH(GZL,GZN,GWF,GWWZ,GZZH,GWWH,
     .             AMSW,GAMW,AMSZ,GAMZ, AMSE,
     .             SPINOR(1,1),SPINOR(1,2),SPINOR(1,3),SPINOR(1,4),
     .             ZVCT,SCH, AMP)
C--
      FACT    = (COLOR*QCD)**(IDP(3,5)-1)
C--
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2)+AMP(3)
     .                  +AMP(4)+AMP(5))**2
      AMP2(3) = FACT*ABS(AMP(6)+AMP(7))**2
C--
C  That's it.
C--
      RETURN
      END
