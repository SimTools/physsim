CC**********************************************************************
C*
C*=================================----===
C* Subroutine FULEETT(IDP,IHEL,X,PV,AMP2)
C*=================================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> e+ e- t tbar.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       X  (1,j) : (R*8) : sin(theta_j/2).
C*          (2,j) : (R*8) : cos(theta_j/2).
C*          (3,j) : (R*8) : phi_j;    (1,2)=(e-,e+).
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  e-.
C*                        :   =  4 for  e+.
C*                        :   =  5 for  t.
C*                        :   =  6 for  tbar.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : non-fusion.
C*        (2)     : (R*8) : fusion.
C*        (3)     : (R*8) : Higgs.
C* (Relation)
C*    Invokes AMPEETT and subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/06/12  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULEETT(IDP,IHEL,X,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      PARAMETER  ( NP = 6 )
      INTEGER*4  IDP(3,NP), IHEL(NP)
      REAL   *8  X(3,2), PV(0:3,NP), AMP2(0:3)
C--
      COMPLEX*16 GTH(2)
      COMPLEX*16 SPINOR(6,NP), AMP(0:3)
      INTEGER*4  NF(NP,2)
      DATA  NF  / +1, -1,   +1, -1,  +1, -1, 
     .            +1, -1,   -1, +1,  -1, +1 /
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
         GWWH      = GW*AMSW
         GZZH      = GW*AMSW/xCOS2W
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
      DO 100 IP = 1, NP
         IF ( NF(IP,2).GT.0 ) THEN
            CALL IXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ELSE
            CALL OXXXXX(PV(0,IP),AMF(IDP(1,IP),IDP(2,IP),IDP(3,IP)),
     .                  IHEL(IP),NF(IP,1), SPINOR(1,IP))
         ENDIF
100   CONTINUE
C--
C  Calculate amplitude**2.
C--
      CALL AMPEETT(GAL,GZL,GZN,GAU,GZU,GWF,GWWA,GWWZ,GW,GZZH,GTH,
     .             AMSW,GAMW,AMSZ,GAMZ,AMSE,AMST,GAMT, AMSH,GAMH,   
     .              	   PV(0,1), PV(0,3),X(1,1),X(2,1),X(3,1),
     .                                           IHEL(1),IHEL(3),
     .                      	    PV(0,4),X(1,2),X(2,2),X(3,2),
     .                                           IHEL(2),IHEL(4),
     .                                 SPINOR(1, 1),SPINOR(1, 2),
     .                                 SPINOR(1, 3),SPINOR(1, 4),
     .                                 SPINOR(1, 5),SPINOR(1, 6),
     .                                                       AMP)
C--
      FACT    = COLOR*QCD
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
