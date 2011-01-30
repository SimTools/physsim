CC**********************************************************************
C*
C*==============================----===
C* Subroutine FULTT(IDP,IHEL,PV,AMP2)
C*==============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar with
C*    decays of t and t-bar into lighter quarks and leptons.
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
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : photon.
C*        (2)     : (R*4) : Z.
C* (Relation)
C*    Invokes AMPTT and subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/07/03  K.Fujii          Original version.
C*    90/08/02  K.Fujii          Modified to accomodate with
C*                               HELAS V1.6.
C*    91/08/01  K.Fujii          Modified to include NR QCD correction.
C*    92/10/27  K.Fujii          Use relativistic off-shell wave
C*                               functions for t and tbar. FQCD is
C*                               now complex. AMPTT is replaced by
C*                               AMPTTP which outputs vector part
C*                               and axial vector part separately.
C*    95/03/06  K.Fujii          Modified to accomodate with 
C*				 HELAS V2.04.
C*    99/08/18  K.Ikematsu       Added P_wave QCD correction.
C*
CC**********************************************************************
 
      SUBROUTINE FULTT(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      INCLUDE 'THRCNS.inc'
#include "SMCONS.inc"
#include "SMPTAB.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2(0:2)
C--
      REAL   *4  GAL(2), GAU(2), GAD(2), GWF(2), GZN(2), GZL(2), GZU(2),
     .           GZD(2), G1(2), GW, GWWZ, GWWA
      REAL   *4  QV(0:3,3)
      COMPLEX*8  SPINOR(6,12), WVCT(0:5,2), AMP(0:2), FQCD, PQCD
      INTEGER*4  NF(8,2)
      DATA  NF  / +1, -1,   -1, -1, +1,   +1, +1, -1,
     .            +1, -1,   +1, +1, -1,   -1, -1, +1 /
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
         AMSH   = AMH
         GAMH   = GMHTOT
         AMSW   = AMW
         GAMW   = GMWTOT
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         CALL UHCUPL(ALF,SN2W, GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1,
     .               GW,GWWZ,GWWA)
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
C  Prepare wave function for W- and W+.
C--
      CALL JIOXXX(SPINOR(1,4),SPINOR(1,5),GWF,AMSW,GAMW,WVCT(0,1))
      CALL JIOXXX(SPINOR(1,8),SPINOR(1,7),GWF,AMSW,GAMW,WVCT(0,2))
C--
C  Prepare spinors for t-bar and t.
C--
      CALL FVIXXX(SPINOR(1,3),WVCT(0,1),GWF,AMST,GAMT,SPINOR(1,11))
      CALL FVOXXX(SPINOR(1,6),WVCT(0,2),GWF,AMST,GAMT,SPINOR(1,12))
C--
C  Calculate QCD correction.
C--
      IF ( NRQCD.EQ.0 ) THEN
         FQCD = 1
         PQCD = 0
      ELSE
         CALL PSUMxx(PV(0,1),PV(0, 2),QV(0,1))
         CALL PSUMxx(PV(0,3),PV(0, 4),QV(0,2))
         CALL PSUMxx(QV(0,2),PV(0, 5),QV(0,2))
         CALL PSUMxx(PV(0,6),PV(0, 7),QV(0,3))
         CALL PSUMxx(QV(0,3),PV(0, 8),QV(0,3))
         CALL RDOTxx(QV(0,3),QV(0,3),S1)
         CALL RDOTxx(QV(0,2),QV(0,2),S2)
         CALL RDOTxx(QV(0,1),QV(0,1),SH)
         RSH = SQRT(SH)
         IF ( IAND(NRQCD,1).NE.0 ) THEN
            CALL QCDMOD(xALFS,AMST,GAMT,AMSH,BETAH,RSH,S1,S2,FQCD)
         ELSE
            FQCD = 1
         ENDIF
         IF ( IAND(NRQCD,2).NE.0 ) THEN
            CALL QCDMOD_P(xALFS,AMST,GAMT,AMSH,BETAH,RSH,S1,S2,PQCD)
         ELSE
            PQCD = 1
         ENDIF
      ENDIF
C--
C  Calculate amplitude**2.
C  AMPTTP outputs
C     AMP(0) = AMP(1) + AMP(2)
C        (1) = v_t part
C        (2) = a_t part.
C  Note that FQCD (S-wave) affects only the v_t part, while PQCD (P-wave)
C  only affects the a_t part.
C--
      CALL AMPTTP(GAL,GAU,GZL,GZU,AMSZ,GAMZ,
     .            SPINOR(1, 1), SPINOR(1, 2),
     .            SPINOR(1,11), SPINOR(1,12), AMP)
      FACT    = 3*COLOR**(IDP(3,5)+IDP(3,7)-2)
     .          * VKM(        3, IDP(1,3), IDP(3,3) )**2
     .          * VKM( IDP(1,4), IDP(1,5), IDP(3,4) )**2
     .          * VKM(        3, IDP(1,6), IDP(3,6) )**2
     .          * VKM( IDP(1,7), IDP(1,8), IDP(3,7) )**2
      AMP2(0) = FACT*ABS(FQCD*AMP(1)+PQCD*AMP(2))**2
      AMP2(1) = FACT*ABS(FQCD*AMP(1))**2
      AMP2(2) = FACT*ABS(PQCD*AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END

