CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULZZF(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> Z + Z with
C*    decays of Z's into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  mb from Z.
C*                        :   =  4 for  m  from Z.
C*                        :   =  5 for  fb from Z.
C*                        :   =  6 for  f  from Z.
C* (Outputs)
C*    AMP2        : (R*4) : (amplitude sum)**2.
C* (Relation)
C*    Invokes AMPZZ and subroutines in HELAS.LOAD.
C* (Update Record)
C*    92/02/24  K.Fujii         Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204. Now uses INSMCN
C*				to initialize /SMCUPL/ to be consistent
C*				with new coupling constant convention.
C*
CC**********************************************************************
 
      SUBROUTINE FULZZF(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2
C--
      REAL   *4  GZF(2,2,2)
      COMPLEX*8  SPINOR(6,6), ZVCT(6,2), AMP(0:2)
      INTEGER*4  NF(6,2)
      DATA  NF  / +1, -1,   -1, +1,   -1, +1,
     .            +1, -1,   +1, -1,   +1, -1   /
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Prepare coupling constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         AMSZ   = AMZ
         GAMZ   = GMZTOT
         AMSE   = AMF(1,2,1)
         GAME   = GMF(1,2,1)
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
      CALL JIOXXX(SPINOR(1,3),SPINOR(1,4),
     .            GZF(1,IDP(2,3),IDP(3,3)),AMSZ,GAMZ,ZVCT(1,1))
C--
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,6),
     .            GZF(1,IDP(2,5),IDP(3,5)),AMSZ,GAMZ,ZVCT(1,2))
C--
C  Calculate amplitude**2.
C--
      CALL AMPZZ (GZL, AMSE,GAME, ZVCT,
     .            SPINOR(1,1), SPINOR(1,2), AMP)
      FACT    = COLOR**(IDP(3,3)+IDP(3,5)-2)
      AMP2    = FACT*ABS(AMP(0))**2
C--
C  That's it.
C--
      RETURN
      END
