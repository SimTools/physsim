CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULTTH(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- --> t + t-bar + g --> t + t-bar + b + b-bar
C*    decays of t, t-bar, b, and b-bar into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  down-type q-bar from t-bar.
C*                        :   =  4 for  up-type f-bar from t-bar.
C*                        :   =  5 for  down-type f from t-bar.
C*                        :   =  6 for  down-type q from t.
C*                        :   =  7 for  up-type f from t.
C*                        :   =  8 for  down-type f-bar from t.
C*                        :   =  9 for  f-bar from gluon.
C*                        :   = 10 for  f from gluon.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C* (Relation)
C*    Invokes ?????? and subroutines in HELAS.LOAD.
C* (Update Record)
C*    2009/06/22  T.Tanabe        First version based on TTHStudy by K.Fujii
C*				
CC**********************************************************************
 
      SUBROUTINE FULTTBB(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,10), IHEL(10)
      REAL   *4  PV(0:3,10), AMP2
C--
      REAL   *4  GFG(2)
C--     
c      REAL   *4  QV(0:3,3)
c      REAL   *4  QVH(0:3)
      COMPLEX*8  SPINOR(6,12), WVCT(6,2), GLUON(6), AMP, FQCD, PQCD
      REAL   *4  CF
      INTEGER*4  NF(10,2)
      DATA  NF  / +1, -1,   -1, -1, +1,   +1, +1, -1,   -1, +1,
     .            +1, -1,   +1, +1, -1,   -1, -1, +1,   +1, -1 /
C--
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
         AMST   = AMF(3,1,2)
         GAMT   = GMF(3,1,2)
         AMSW   = AMW
         GAMW   = GMWTOT
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         AMSH   = AMH
         GAMH   = GMHTOT
C--
c        set QCD coupling for ttg
         GFG(1) =  -SQRT(4d0*xPI*xALFS)
         GFG(2) =   GFG(1)
C--
         QCD    = (1+xALFS/xPI)
         COLOR  = 3*QCD
         CF     = 2d0
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
C  Prepare wave function for W- and W+.
C--
      CALL JIOXXX(SPINOR(1,4),SPINOR(1,5),GWF,AMSW,GAMW,WVCT(1,1))
      CALL JIOXXX(SPINOR(1,8),SPINOR(1,7),GWF,AMSW,GAMW,WVCT(1,2))
C--
C  Prepare wave function for gluon.
C--
      CALL JIOXXX(SPINOR(1,9),SPINOR(1,10),GFG,ZERO,ZERO,GLUON)
C--
C  Prepare spinors for t-bar and t.
C--
      CALL FVIXXX(SPINOR(1,3),WVCT(1,1),GWF,AMST,GAMT,SPINOR(1,11))
      CALL FVOXXX(SPINOR(1,6),WVCT(1,2),GWF,AMST,GAMT,SPINOR(1,12))
C--
C  Calculate QCD correction.
C--
c#ifdef TT_THRESHOLD
cC      WRITE(*,*)'correction'
c      CALL PSUMxx(PV(0,4),PV(0, 5),QV(0,1))
c      CALL PSUMxx(QV(0,1),PV(0, 3),QV(0,1))
c      CALL PSUMxx(PV(0,7),PV(0, 8),QV(0,2))
c      CALL PSUMxx(QV(0,2),PV(0, 6),QV(0,2))
c      CALL PSUMxx(QV(0,1),QV(0, 2),QV(0,3))
c      CALL RDOTxx(QV(0,1),QV(0,1),S1)
c      CALL RDOTxx(QV(0,2),QV(0,2),S2)
c      CALL RDOTxx(QV(0,3),QV(0,3),SH)
c      RSH = SQRT(SH)
c      CALL QCDMOD(xALFS,AMST,GAMT,AMSH,BETAH,RSH,S1,S2,FQCD)
cC--
cC
cC   no P-wave correction 
cC
cC--
c      PQCD = 1
cC     CALL QCDMOD_P(xALFS,AMST,GAMT,AMSH,BETAH,RSH,S1,S2,PQCD)
c#else
cC      WRITE(*,*)'no correction'      
      FQCD = 1
      PQCD = 1
cC         WRITE(*,*)'no correction'
c#endif
C--
C  Calculate amplitude**2.
C--
      CALL AMPTTG(GAL,GAU,GZL,GZU,GFG, AMSZ,GAMZ,AMST,GAMT,
     .            SPINOR(1, 1), SPINOR(1, 2),
     .            SPINOR(1,11), SPINOR(1,12), GLUON, AMP)
C--
C  No QCD correction to glue --> ffbar at present.
C--
c      FACT    = COLOR**(IDP(3,5)+IDP(3,7)-1)
c     .          * 3**(IDP(3,9)-1)
      FACT    = COLOR**(IDP(3,5)+IDP(3,7)-2)
     .          * CF
     .          * VKM(        3, IDP(1,3), IDP(3,3) )**2
     .          * VKM( IDP(1,4), IDP(1,5), IDP(3,4) )**2
     .          * VKM(        3, IDP(1,6), IDP(3,6) )**2
     .          * VKM( IDP(1,7), IDP(1,8), IDP(3,7) )**2
      AMP2 = FACT*ABS(FQCD*AMP)**2
C--
C  That's it.
C--
      RETURN
      END
