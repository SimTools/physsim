CC**********************************************************************
C*
C*===============================----===
C* Subroutine FULSEE(IDP,IHEL,PV,AMP2)
C*===============================----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> SE+ + SE-.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*       PV(*, i) : (R*4) : 4-momentum of i-th particle.
C*                        : i =  1 for  e-.
C*                        :   =  2 for  e+.
C*                        :   =  3 for  LSP  from SE-.
C*                        :   =  4 for  ld   from SE-.
C*                        :   =  5 for  LSP  from SE+.
C*                        :   =  6 for  ldb  from SE+.
C* (Outputs)
C*    AMP2(0)     : (R*4) : (amplitude sum)**2.
C*        (1)     : (R*4) : photon.
C*        (2)     : (R*4) : Z.
C*        (3)     : (R*4) : t-channel X0-exchange.
C* (Relation)
C*    Invokes AMPTT and subroutines in HELAS.LOAD.
C* (Update Record)
C*    92/06/28  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE FULSEE(IDP,IHEL,PV,AMP2)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
#include "SFSFCN.inc"
#include "SMCONS.inc"
#include "SSCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SMCUPL.inc"
#include "SSCUPL.inc"
      INTEGER*4  IDP(3,*), IHEL(*)
      REAL   *4  PV(0:3,*), AMP2(0:4)
C--
      COMPLEX*8  SPINOR(6,20), EIN(6), EOT(6), SFIN(3), SFOT(3),
     .           AMP(0:6)
      COMPLEX*8  CX(2,2,4), GNESE(2,2,4), GNESEC(2,2,4)
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
         AMSZ      = AMZ
         GAMZ      = GMZTOT
C--
C  Select final sate handedness.
C     IE/IP = 1 : se_L
C           = 2 : se_R
C  where IE and IP correspond to selectron and spositron.
C--
         IE        = IHANDM
         IP        = IHANDP
C--
C  Initialize coupling constants.
C--
         IF ( IE.EQ.IP ) THEN
            GAX    = GAL(1)
            GZX    = GZL(IE)
         ELSE
            GAX    = 0
            GZX    = 0
         ENDIF
C--
         DO 10 IX = 1, 4
            CX(1,1,IX) = GNESEL(1,IX)
            CX(2,1,IX) = GNESEL(2,IX)
            CX(1,2,IX) = GNESER(1,IX)
            CX(2,2,IX) = GNESER(2,IX)
10       CONTINUE
         DO 20 K = 1, 4
            GNESE(1,1,K)  = CX(1,IE,K)
            GNESE(2,1,K)  = CX(2,IE,K)
            GNESE(1,2,K)  = CONJG(CX(2,IP,K))
            GNESE(2,2,K)  = CONJG(CX(1,IP,K))
            GNESEC(1,1,K) = CONJG(CX(2,IE,K))
            GNESEC(2,1,K) = CONJG(CX(1,IE,K))
            GNESEC(1,2,K) = CX(1,IP,K)
            GNESEC(2,2,K) = CX(2,IP,K)
20       CONTINUE
      ENDIF
C--
C  Prepare spinors for external lines.
C--
      CALL IXXXXX(PV(0,1),AMSE  ,IHEL(1),+1,EIN)
      CALL OXXXXX(PV(0,2),AMSE  ,IHEL(2),-1,EOT)
C--
      AM4 = AMF(IDP(1,4),IDP(2,4),IDP(3,4))
      AM6 = AMF(IDP(1,6),IDP(2,6),IDP(3,6))
      CALL IXXXXX(PV(0,3),SZM(1),IHEL(3),-1,SPINOR(1,3))
      CALL OXXXXX(PV(0,4),AM4   ,IHEL(4),+1,SPINOR(1,4))
      CALL OXXXXX(PV(0,5),SZM(1),IHEL(5),+1,SPINOR(1,5))
      CALL IXXXXX(PV(0,6),AM6   ,IHEL(6),-1,SPINOR(1,6))
C--
C  Prepare spinors for S+ and S-.
C--
      CALL HIOXXX(SPINOR(1,3),SPINOR(1,4),GNESEC(1,1,1),
     .            SFM(IE+1),GMSF(IE+1), SFOT )
      CALL HIOXXX(SPINOR(1,6),SPINOR(1,5),GNESEC(1,2,1),
     .            SFM(IP+1),GMSF(IP+1), SFIN )
C--
C  Calculate amplitude**2.
C--
      CALL AMPSSE(GAL,GAX,GZL,GZX,AMSZ,GAMZ,SZM,GMSZ,GNESE,
     .                   EIN,EOT,SFIN,SFOT,AMP)
      AMP2(0) = ABS(AMP(0))**2
      AMP2(1) = ABS(AMP(1))**2
      AMP2(2) = ABS(AMP(2))**2
      AMP2(3) = ABS(AMP(3)+AMP(4)+AMP(5)+AMP(6))**2
      AMP2(4) = ABS(AMP(1)+AMP(2))**2
C--
C  That's it.
C--
      RETURN
      END
