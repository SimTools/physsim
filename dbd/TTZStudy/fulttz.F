CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULTTZ(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar + Z with
C*    decays of t, t-bar, and Z into lighter quarks and leptons.
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
C*                        :   =  9 for  f-bar from Z.
C*                        :   = 10 for  f from Z.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : amp**2 of non-Higgs diagrams.
C*        (2)     : (R*4) : amp**2 of Higgs diagram.
C* (Relation)
C*    Invokes AMPTTZ and subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/06/15  K.Fujii         Original version.
C*    90/08/02  K.Fujii         Modified to accomodate with V1.06.
C*    95/03/16  K.Fujii         Modified to accomodate with V2.04.
C*				Notice that GTH(2) is now complex and
C*                              GZH has an opposite sign.
C*				UHCUPL must be the one in xcxc_lib.
C*    95/06/09  K.Fujii         Removed UHCUPL call. Now use /SMCUPL/.
C*				Added QCD correction to ttZ.
C*    09/08/13  K.Ikematsu      Implemented S-wave ttbar threshold
C*                              enhancement
C*
CC**********************************************************************
 
      SUBROUTINE FULTTZ(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2(0:2)
C--
      REAL   *4  GZH
      COMPLEX*8  GTH(2)
C--
      REAL   *4  QV(0:3,3)
      COMPLEX*8  FQCD
C--     
      REAL   *4  GZF(2,2,2)
      COMPLEX*8  SPINOR(6,12), WVCT(0:5,2), ZVCT(0:5), AMP(0:2)
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
         AMSW   = AMW
         GAMW   = GMWTOT
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         AMSH   = AMH
         GAMH   = GMHTOT
C--
         GZH    =   GW*AMSW/CS2W
C--
         GZF(1,1,1) = GZN(1)
         GZF(2,1,1) = GZN(2)
         GZF(1,2,1) = GZL(1)
         GZF(2,2,1) = GZL(2)
         GZF(1,1,2) = GZU(1)
         GZF(2,1,2) = GZU(2)
         GZF(1,2,2) = GZD(1)
         GZF(2,2,2) = GZD(2)
C--
         COLOR      = 3*(1+xALFS/xPI)
      ENDIF
C--
      AMST   =   AMF(3,1,2)
      GAMT   =   GMF(3,1,2)
      GTH(1) = - GW*AMST/AMSW/2
      GTH(2) =   GTH(1)
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
      CALL JIOXXX(SPINOR(1,4),SPINOR(1,5),GWF,AMSW,GAMW,WVCT(0,1))
      CALL JIOXXX(SPINOR(1,8),SPINOR(1,7),GWF,AMSW,GAMW,WVCT(0,2))
C--
C  Prepare wave function for Z.
C--
      CALL JIOXXX(SPINOR(1,9),SPINOR(1,10),
     .            GZF(1,IDP(2,9),IDP(3,9)), AMSZ,GAMZ, ZVCT)
C--
C  Prepare spinors for t-bar and t.
C--
      CALL FVIXXX(SPINOR(1,3),WVCT(0,1),GWF,AMST,GAMT,SPINOR(1,11))
      CALL FVOXXX(SPINOR(1,6),WVCT(0,2),GWF,AMST,GAMT,SPINOR(1,12))
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
#else
C      WRITE(*,*)'no correction'      
      FQCD = 1
#endif
C--
C  Calculate amplitude**2.
C--
      CALL AMPTTZ(GAL,GAU,GZL,GZU,GZH,GTH,
     .            AMSZ,GAMZ,AMSE,AMST,GAMT,AMSH,GAMH, ZVCT,
     .            SPINOR(1, 1), SPINOR(1, 2),
     .            SPINOR(1,11), SPINOR(1,12), AMP)
      FACT    = COLOR**(IDP(3,5)+IDP(3,7)+IDP(3,9)-2)
     .          * VKM(        3, IDP(1,3), IDP(3,3) )**2
     .          * VKM( IDP(1,4), IDP(1,5), IDP(3,4) )**2
     .          * VKM(        3, IDP(1,6), IDP(3,6) )**2
     .          * VKM( IDP(1,7), IDP(1,8), IDP(3,7) )**2
      AMP2(0) = FACT*ABS(FQCD*AMP(1)+AMP(2))**2
      AMP2(1) = FACT*ABS(FQCD*AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END
