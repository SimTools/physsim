CC**********************************************************************
C*
C*==============================----===
C* Subroutine FULWW(IDP,IHEL,PV,AMP2)
C*==============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> W- + W+  with
C*    decays of W- and  W+ Z into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  down-type f from W-.
C*                        :   =  4 for  up-type f-bar from W-.
C*                        :   =  5 for  down-type f-bar from W+.
C*                        :   =  6 for  up-type f from W+.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : amp**2 of s-channel diagrams.
C*        (2)     : (R*4) : amp**2 of t-channel diagram.
C* (Relation)
C*    Invokes AMPWW and subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/08/13  K.Fujii         Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204. Now uses INSMCN
C*				to initialize /SMCUPL/ to be consistent
C*				with new coupling constant convention.
C*
CC**********************************************************************
 
      SUBROUTINE FULWW(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,6), IHEL(6)
      REAL   *4  PV(0:3,6), AMP2(0:2)
C--
      COMPLEX*8  SPINOR(6,6), WVCT(6,2), AMP(0:2)
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
         NCALL  = 1
         AMSW   = AMW
         GAMW   = GMWTOT
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         COLOR  = 3*(1+xALFS/xPI)
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
C  Prepare wave function for W- and W+.
C--
      CALL JIOXXX(SPINOR(1,4),SPINOR(1,3),GWF,AMSW,GAMW,WVCT(1,1))
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,6),GWF,AMSW,GAMW,WVCT(1,2))
C--
C  Calculate amplitude**2.
C--
      CALL AMPWW(GAL,GZL,GWF,GW,AMSZ,GAMZ,
     .           WVCT,SPINOR(1,1),SPINOR(1,2), AMP)
      FACT    = COLOR**(IDP(3,3)+IDP(3,5)-2)
     .          * VKM( IDP(1,4), IDP(1,3), IDP(3,3) )**2
     .          * VKM( IDP(1,6), IDP(1,5), IDP(3,5) )**2
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END
