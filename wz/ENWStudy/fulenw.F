CC**********************************************************************
C*
C*=================================----===
C* Subroutine FULENW(IDP,IHEL,X,PV,AMP2)
C*=================================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> e- + nubar + W+.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       X   ( 1) : (R*8) : sin(theta_e/2).
C*           ( 2) : (R*8) : cos(theta_e/2).
C*           ( 3) : (R*8) : phi_e.
C*       PV(*, i) : (R*8) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  e-.
C*                        :   =  4 for  nubar.
C*                        :   =  5 for  fdb from W+.
C*                        :   =  6 for  fu from W+.
C* (Outputs)
C*    AMP2(0)     : (R*8) : (amplitude sum)**2.
C*        (1)     : (R*8) : photon.
C*        (2)     : (R*8) : Z.
C*        (3)     : (R*8) : W.
C* (Relation)
C*    Invokes AMPENW and subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/04/22  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULENW(IDP,IHEL,X,PV,AMP2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SMCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *8  X(3), PV(0:3,*), AMP2(0:3)
C--
      COMPLEX*16 SPINOR(6,6), WVCT(6), AMP(0:6)
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
         NCALL     = 1
         AMSE      = AMF(1,2,1)
         AMSW      = AMW
         GAMW      = GMWTOT
         AMSZ      = AMZ
         GAMZ      = GMZTOT
C--
         COLOR     = 3
         QCD       = 1 + xALFS/xPI
C--
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
C  Prepare wave function for W+.
C--
      CALL JIOXXX(SPINOR(1,5),SPINOR(1,6),GWF,AMSW,GAMW,WVCT(1))
C--
C  Calculate amplitude**2.
C--
      CALL AMPENW(GAL,GZL,GZN,GWF,GWWA,GWWZ,AMSZ,GAMZ,AMSW,GAMW,AMSE,
     .		  PV(0,1),PV(0,3),X(1),X(2),X(3), IHEL(1),IHEL(3),
     .            SPINOR(1,1),SPINOR(1,2),SPINOR(1,3),SPINOR(1,4),
     .		  WVCT(1),AMP)
C--
      FACT    = (COLOR*QCD)**(IDP(3,5)-1)
     .          * VKM( IDP(1,6), IDP(1,5), IDP(3,5) )**2
C--
      AMP2(0) = FACT*ABS(AMP(0)              )**2
      AMP2(1) = FACT*ABS(AMP(1)+AMP(3)       )**2
      AMP2(2) = FACT*ABS(AMP(2)+AMP(4)+AMP(5))**2
      AMP2(3) = FACT*ABS(AMP(6)              )**2
C--
C  That's it.
C--
      RETURN
      END
