CC**********************************************************************
C*
C*=================================----===
C* Subroutine FULNNZZ(IDP,IHEL,X,PV,AMP2)
C*=================================----===
C*
C* (Purpose)
C*    Calculate amplitude for e- e+ ---> nu nub Z Z.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       X  (1,j) : (R*8) : sin(theta_j/2).
C*          (2,j) : (R*8) : cos(theta_j/2).
C*          (3,j) : (R*8) : phi_j;    (1,2)=(nu,nub).
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  nu.
C*                        :   =  4 for  nub.
C*                        :   =  5 for  f  from Z1.
C*                        :   =  6 for  fb from Z1.
C*                        :   =  7 for  f  from Z2.
C*                        :   =  8 for  fb from Z2.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : non-fusion.
C*        (2)     : (R*8) : fusion.
C*        (3)     : (R*8) : Higgs.
C* (Relation)
C*    Invokes AMPNNZZ and subroutines in HELAS.LOAD.
C* (Update Record)
C*    2005/03/06  K.Fujii        Original JSF version.
C*
CC**********************************************************************
 
      SUBROUTINE FULNNZZ(IDP,IHEL,X,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *8  X(3,2), PV(0:3,*), AMP2(0:3)
C--
      REAL   *8  GZF(2,2,2)
      COMPLEX*16 SPINOR(6,8), VCT(6,2), AMP(0:3)
      INTEGER*4  NF(8,2)
      DATA  NF  / +1, -1,   +1, -1,   +1, -1, +1, -1,
     .            +1, -1,   -1, +1,   -1, +1, -1, +1 /
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
         AMSW      = AMW
         GAMW      = GMWTOT
         AMSZ      = AMZ
         GAMZ      = GMZTOT
         AMSH      = AMH
         GAMH      = GMHTOT
C--
         GWWH      = GW*AMSW
         GZZH      = GW*AMSW/xCOS2W
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
         COLOR     = 3
         QCD       = 1 + xALFS/xPI
C--
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
C  Prepare wave function for Z's.
C--
      CALL JIOXXX(SPINOR(1,6),SPINOR(1,5),
     .            GZF(1,IDP(2,5),IDP(3,5)),AMSZ,GAMZ,VCT(1,1))
      CALL JIOXXX(SPINOR(1,8),SPINOR(1,7),
     .            GZF(1,IDP(2,7),IDP(3,7)),AMSZ,GAMZ,VCT(1,2))
C--
C  Calculate amplitude**2.
C--
      CALL AMPNNZZ(GZL,GZN,GWF,GWWZ,GWWH,GZZH,
     .              AMSW,GAMW, AMSZ,GAMZ, AMSE, AMSH,GAMH,   
     .                            SPINOR(1,1),SPINOR(1,2),
     .                            SPINOR(1,3),SPINOR(1,4),
     .                                           VCT, AMP)
C--
      FACT    = (COLOR*QCD)**(IDP(3,5)+IDP(3,7)-2)
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
