CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULNNTT(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> nu nubar t tbar.
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
C*                        :   =  5 for  b   from t.
C*                        :   =  6 for  fdb from W+.
C*                        :   =  7 for  fu  from W+.
C*                        :   =  8 for  bb  from t.
C*                        :   =  9 for  fd  from W-.
C*                        :   = 10 for  fub from W-.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : non-fusion.
C*        (2)     : (R*8) : fusion.
C*        (3)     : (R*8) : Higgs.
C* (Relation)
C*    Invokes AMPNNTT and subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/06/09  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULNNTT(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *8  PV(0:3,*), AMP2(0:3)
C--
      COMPLEX*16 GTH(2)
      COMPLEX*16 SPINOR(6,12), VCT(6,2), AMP(0:3)
      INTEGER*4  NF(10,2)
      DATA  NF  / +1, -1,   +1, -1,  +1, -1, +1,  -1, +1, -1,
     .            +1, -1,   -1, +1,  -1, +1, -1,  +1, -1, +1 /
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
         AMST      = AMF(3,1,2)
         GAMT      = GMF(3,1,2)
         AMSB      = AMF(3,2,2)
         GAMB      = GMF(3,2,2)
         AMSW      = AMW
         GAMW      = GMWTOT
         AMSZ      = AMZ
         GAMZ      = GMZTOT
         AMSH      = AMH
         GAMH      = GMHTOT
C--
         GWWH      =   GW*AMSW
         GTH(1)    = - GW*AMST/AMSW/2
         GTH(2)    =   GTH(1)
C--
         COLOR     = 3
         QCD       = 1 + xALFS/xPI
C--
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      DO 100 IP = 1, 10
         IF ( NF(IP,2).GT.0 ) THEN
            CALL IXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ELSE
            CALL OXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ENDIF
100   CONTINUE
C--
C  Prepare wave function for W's.
C--
      CALL JIOXXX(SPINOR(1, 6),SPINOR(1,7),GWF,AMSW,GAMW,VCT(1,1))
      CALL JIOXXX(SPINOR(1,10),SPINOR(1,9),GWF,AMSW,GAMW,VCT(1,2))
C--
C  Prepare spinors for t-bar and t.
C--
      CALL FVOXXX(SPINOR(1,5),VCT(1,1),GWF,AMST,GAMT,SPINOR(1,11))
      CALL FVIXXX(SPINOR(1,8),VCT(1,2),GWF,AMST,GAMT,SPINOR(1,12))
C--
C  Calculate amplitude**2.
C--
      CALL AMPNNTT(GAL,GZL,GZN,GAU,GZU,GWF, GW, GWWH, GTH,
     .              AMSW,GAMW, AMSZ,GAMZ, AMSE, AMSB,GAMB,
     .              AMSH,GAMH,  SPINOR(1, 1),SPINOR(1, 2),
     .                          SPINOR(1, 3),SPINOR(1, 4),
     .                          SPINOR(1,11),SPINOR(1,12), AMP)
C--
      FACT    = (COLOR*QCD)**(IDP(3,6)+IDP(3,9)-1)
     .          * VKM(         3, IDP(1,5), IDP(3,5) )**2
     .          * VKM( IDP(1, 7), IDP(1,6), IDP(3,6) )**2
     .          * VKM(         3, IDP(1,8), IDP(3,8) )**2
     .          * VKM( IDP(1,10), IDP(1,9), IDP(3,9) )**2
C--
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
      AMP2(3) = FACT*ABS(AMP(3))**2
C--
C  That's it.
C--
      RETURN
      END
