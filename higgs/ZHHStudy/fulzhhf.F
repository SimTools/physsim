CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULZHHF(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> Z + H + H with
C*    decays of Z and H into lighter quarks and leptons.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  fb from Z.
C*                        :   =  4 for  f  from Z.
C*                        :   =  5 for  H1.
C*                        :   =  6 for  H2.
C* (Outputs)
C*    AMP2        : (R*4) : (amplitude sum)**2.
C* (Relation)
C*    Invokes AMPZHH and subroutines in HELAS.LOAD.
C* (Update Record)
C*    08/11/30  K.Fujii         Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULZHHF(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2
C--
      REAL   *4  GZZH, GHHH, GZZHH
C--
      REAL   *4  QV(0:3,6), GZF(2,2,2)
C--
      COMPLEX*8  SPINOR(6,6), ZVCT(6), HIGG1(3), HIGG2(3), AMP(0:3)
      INTEGER*4  NF(4,2)
      DATA  NF  / +1, -1,   -1, +1,
     .            +1, -1,   +1, -1 /
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
         AMSH   = AMH
         GAMH   = GMHTOT
         AMSW   = AMW
         AMSZ   = AMZ
         GAMZ   = GMZTOT
C--
         V       =   2*AMSW/GW
         GZZH    =   GW*AMSW/CS2W
         GHHH    =  -AMSH**2/V * 3 
         GZZHH   =   GW**2/xCOS2W/2
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
      DO 100 IP = 1, 4
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
      CALL SXXXXX(PV(0,5),+1,HIGG1)
      CALL SXXXXX(PV(0,6),+1,HIGG2)
C--
C  Prepare wave function for Z.
C--
      CALL JIOXXX(SPINOR(1,3),SPINOR(1,4),
     .            GZF(1,IDP(2,3),IDP(3,3)), AMSZ,GAMZ, ZVCT)
C--
C  Calculate amplitude**2.
C--
      CALL AMPZHH(GZZH,GHHH,GZZHH,GZL, AMSZ,GAMZ, AMSH, GAMH,
     .            ZVCT,HIGG1,HIGG2, 
     .            SPINOR(1,1),SPINOR(1,2), AMP)
      FACT    = (COLOR*QCD)**(IDP(3,3)-1)
      AMP2    = FACT*ABS(AMP(0))**2
C--
C  That's it.
C--
      RETURN
      END
