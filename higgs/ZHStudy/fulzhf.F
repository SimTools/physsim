CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULZHF(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> Z + H with
C*    decays of Z and H into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  bb from H.
C*                        :   =  4 for  b  from H.
C*                        :   =  5 for  mb from Z.
C*                        :   =  6 for  m  from Z.
C* (Outputs)
C*    AMP2        : (R*4) : (amplitude sum)**2.
C* (Relation)
C*    Invokes AMPZH and subroutines in HELAS.LOAD.
C* (Update Record)
C*    92/02/24  K.Fujii         Original version.
C*    95/03/25  K.Fujii         Use INSMCN instead of CONINT.
C*				Notice that GBH(2) is now complex and
C*                              GZH has an opposite sign.
C*				UHCUPL must be the one in xcxc_lib.
C*
CC**********************************************************************
 
      SUBROUTINE FULZHF(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2
C--
      REAL   *4  GZH
      COMPLEX*8  GBH(2)
C--
      REAL   *4  QV(0:3,6), GZF(2,2,2)
C--
      COMPLEX*8  SPINOR(6,6), ZVCT(6), HIGG(3), AMP
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
         CS2W   = xCOS2W
         AMSB   = AMF(3,2,2)
         AMSH   = AMH
         GAMH   = GMHTOT
         AMSW   = AMW
         AMSZ   = AMZ
         GAMZ   = GMZTOT
C--
         GZH    =   GW*AMSW/CS2W
         GBH(1) = - GW*AMSB/AMSW/2
         GBH(2) = GBH(1)
         GZF(1,1,1) = GZN(1)
         GZF(2,1,1) = GZN(2)
         GZF(1,2,1) = GZL(1)
         GZF(2,2,1) = GZL(2)
         GZF(1,1,2) = GZU(1)
         GZF(2,1,2) = GZU(2)
         GZF(1,2,2) = GZD(1)
         GZF(2,2,2) = GZD(2)
         COLOR      = 3
         QCD        = 1 + xALFS/xPI
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
C  Prepare wave function for H.
C--
      CALL HIOXXX(SPINOR(1,3),SPINOR(1,4),GBH,AMSH,GAMH,HIGG)
C--
C  Prepare wave function for Z.
C--
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,6),
     .            GZF(1,IDP(2,5),IDP(3,5)), AMSZ,GAMZ, ZVCT)
C--
C  Calculate amplitude**2.
C--
      CALL AMPZH (GZH,GZL, AMSZ,GAMZ, ZVCT,HIGG, 
     .            SPINOR(1,1),SPINOR(1,2), AMP)
      FACT    = COLOR**(IDP(3,3)-1)*(COLOR*QCD)**(IDP(3,5)-1)
      AMP2    = FACT*ABS(AMP)**2
C--
C  That's it.
C--
      RETURN
      END
