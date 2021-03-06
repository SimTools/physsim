CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULTTZ(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar + Z.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  f-bar.
C*                        :   =  4 for  f.
C*                        :   =  5 for  Z.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : amp**2 of non-Higgs diagrams.
C*        (2)     : (R*4) : amp**2 of Higgs diagram.
C* (Relation)
C*    Invokes AMPTTZ and subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/06/09  K.Fujii         Original version.
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
      COMPLEX*8  GFH(2)
C--     
      COMPLEX*8  SPINOR(6,4), ZVCT(0:5), AMP(0:2)
      INTEGER*4  NF(4,2)
      DATA  NF  / +1, -1,   -1, +1, 
     .            +1, -1,   +1, -1 /
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
         COLOR      = 3*(1+xALFS/xPI)
      ENDIF
C--
      AMSF   =   AMF(IDP(1,3),IDP(2,3),IDP(3,3))
      GAMF   =   GMF(IDP(1,3),IDP(2,3),IDP(3,3))
      GFH(1) = - GW*AMSF/AMSW/2
      GFH(2) =   GFH(1)
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
C  Prepare wave function for Z.
C--
      CALL VXXXXX(PV(0,5),AMSZ,IHEL(5),+1,ZVCT)
C--
C  Calculate amplitude**2.
C--
      CALL AMPTTZ(GAL,GAU,GZL,GZU,GZH,GFH,
     .            AMSZ,GAMZ,AMSE,AMSF,GAMF,AMSH,GAMH, ZVCT,
     .            SPINOR(1,1), SPINOR(1,2),
     .            SPINOR(1,3), SPINOR(1,4), AMP)
      FACT    = COLOR**(IDP(3,3)-1)
      AMP2(0) = FACT*ABS(AMP(0))**2
      AMP2(1) = FACT*ABS(AMP(1))**2
      AMP2(2) = FACT*ABS(AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END
