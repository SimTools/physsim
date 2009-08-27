CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULTTH(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar + H_SM with
C*    decays of t, t-bar, and H_SM into lighter quarks and leptons.
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
C*                        :   =  9 for  f-bar from H_SM.
C*                        :   = 10 for  f from H_SM.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : amp**2 of Higgsostrahlung diagram.
C*        (2)     : (R*4) : amp**2 of ZH diagram.
C* (Relation)
C*    Invokes AMPTTH and subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/03/28  K.Fujii         Original version.
C*    07/11/07  K.Fujii         Allow no Higgs decay
C*    07/11/08  A.Ishikawa      Allow ttbar threshold enhancement
C*				
CC**********************************************************************
 
      SUBROUTINE FULTTH(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,10), IHEL(10)
      REAL   *4  PV(0:3,10), AMP2(0:2)
C--
      REAL   *4  GZH
      COMPLEX*8  GTH(2), GFH(2)
C--     
      REAL   *4  QV(0:3,3)
      REAL   *4  QVH(0:3)
      COMPLEX*8  SPINOR(6,12), WVCT(6,2), HIGG(3), AMP(0:2), FQCD, PQCD
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
         GZH    =   GW*AMSW/CS2W
         GTH(1) = - GW*AMST/AMSW/2
         GTH(2) =   GTH(1)
C--
         QCD    = (1+xALFS/xPI)
         COLOR  = 3*QCD
      ENDIF
C--
      GFH(1) = - GW*AMF(IDP(1,9),IDP(2,9),IDP(3,9))/AMSW/2
      GFH(2) =   GFH(1)
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
C  Prepare wave function for H_SM.
C--
#ifdef NO_H_DECAY
      CALL PSUMXX(PV(0,9),PV(0,10),QVH)
      CALL SXXXXX(QVH,1,HIGG)
#else
      CALL HIOXXX(SPINOR(1,9),SPINOR(1,10),GFH,AMSH,GAMH,HIGG)
#endif
C--
C  Prepare spinors for t-bar and t.
C--
      CALL FVIXXX(SPINOR(1,3),WVCT(1,1),GWF,AMST,GAMT,SPINOR(1,11))
      CALL FVOXXX(SPINOR(1,6),WVCT(1,2),GWF,AMST,GAMT,SPINOR(1,12))
C--
C  Calculate QCD correction.
C--
#ifdef TT_THRESHOLD
C      WRITE(*,*)'correction'
      CALL PSUMxx(PV(0,4),PV(0, 5),QV(0,1))
      CALL PSUMxx(QV(0,1),PV(0, 3),QV(0,1))
      CALL PSUMxx(PV(0,7),PV(0, 8),QV(0,2))
      CALL PSUMxx(QV(0,2),PV(0, 6),QV(0,2))
      CALL PSUMxx(QV(0,1),QV(0, 2),QV(0,3))
      CALL RDOTxx(QV(0,1),QV(0,1),S1)
      CALL RDOTxx(QV(0,2),QV(0,2),S2)
      CALL RDOTxx(QV(0,3),QV(0,3),SH)
      RSH = SQRT(SH)
      CALL QCDMOD(xALFS,AMST,GAMT,AMSH,BETAH,RSH,S1,S2,FQCD)
C--
C
C   no P-wave correction 
C
C--
      PQCD = 1
C     CALL QCDMOD_P(xALFS,AMST,GAMT,AMSH,BETAH,RSH,S1,S2,PQCD)
#else
C      WRITE(*,*)'no correction'      
      FQCD = 1
      PQCD = 1
C         WRITE(*,*)'no correction'
#endif
C--
C  Calculate amplitude**2.
C--
      CALL AMPTTH(GAL,GAU,GZL,GZU,GZH,GTH, AMSZ,GAMZ,AMST,GAMT,
     .            SPINOR(1, 1), SPINOR(1, 2),
     .            SPINOR(1,11), SPINOR(1,12), HIGG, AMP)
C--
C  No QCD correction to H --> ffbar at present.
C--
      FACT    = COLOR**(IDP(3,5)+IDP(3,7)-1)
     .          * 3**(IDP(3,9)-1)
     .          * VKM(        3, IDP(1,3), IDP(3,3) )**2
     .          * VKM( IDP(1,4), IDP(1,5), IDP(3,4) )**2
     .          * VKM(        3, IDP(1,6), IDP(3,6) )**2
     .          * VKM( IDP(1,7), IDP(1,8), IDP(3,7) )**2
      AMP2(0) = FACT*ABS(FQCD*AMP(1)+FQCD*AMP(2))**2
      AMP2(1) = FACT*ABS(FQCD*AMP(1))**2
      AMP2(2) = FACT*ABS(FQCD*AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END
