C ==============================================================
      SUBROUTINE ISR_WHIZARD_SPECTRUM(Z1, Z2, ZE1, ZE2, ZF1, ZF2,
     >  xALF0, LORDER, ROOTS, EMBM, EPBM, QV, PISR, WAT, RS) 
C
C(Function) 
C  Generate ISR spectrum and photon four momenta a la Whizard.
C(Input)
C  Z1,Z2: random varaiable for ISR spectrum
C  ZE1,ZE2: random variable to generate gamma Pt
C  ZF1,ZF2: random variable to generate gamma azimuthal angle
C  xALF0  : fine structure constant (1/137....)
C  LORDER : LLA_Order, 1, 2, or 3.
C  ROOTS  : Nominal CM energy
C  EMBM   : e- energy after beam strahlung
C  EPBM   : e+ energy after beam strahlung.
C  WAT    : initial weight of this sample
C(Output)
C  QV(0:3) : CMS four vector after ISR
C  PISR(0:3,2) : four vector of 2 ISR gamma
C  WAT     : WAT of this sample.
C  RS      :  New center of mass energy, < 0 if too low.
C
C(Author)
C  Akiya Miyamoto, 21-June-2011
C
C ---------------------------------------------------------------
      IMPLICIT NONE
      REAL*8   Z1, Z2, ZE1, ZE2, ZF1, ZF2
      REAL*8   PISR(0:3,2), WAT
      REAL*4   xALF0, EMBM, EPBM, QV(0:3) 
      REAL*8   EPSISR, FISR3, FISR4, Z3, Z4, Z3INI, Z4INI, AME
      INTEGER*4  LORDER,KK
      REAL*4   ROOTS, RS, RSSQ, PA
#if 0
#include "SMCONS.inc"
#include "SMPTAB.inc"

C  Generate ISR spectrum a la Whizard
       AME = AMF(1,2,1)
#else
       REAL*8  xPI
       INTEGER*4 NCALL
       DATA NCALL /0/
       IF (NCALL .EQ.0) THEN
          NCALL = 1
          xPI = DACOS(-1.D0)
          AME = 0.5110034D-3
       ENDIF
#endif
          EPSISR=xALF0/xPI*2*log(RS/AME)

          Z3INI=EMBM/ROOTS*2
          Z3   =Z1
          FISR3=1.0D0
          CALL ISR_FUNCTION(FISR3, Z3, EPSISR, LORDER)
          Z4INI=EPBM/ROOTS*2
          Z4   =Z2
          FISR4=1.0D0
          CALL ISR_FUNCTION(FISR4, Z4, EPSISR, LORDER)
          WAT=WAT*FISR3*FISR4

         CALL ISR_REMNANT(Z3,Z3INI,ZE1,ZF1,DBLE(ROOTS),PISR(0,1))
         CALL ISR_REMNANT(Z4,Z4INI,ZE2,ZF2,DBLE(ROOTS),PISR(0,2))
C
C PISRDBL: 4 vector of ISR photons
C QV(i,1): 4 vector of initial state. Recalculate RS for ME calculation 
C   
          DO 10510 KK=1,3
            PISR(KK,2)=-PISR(KK,2)
10510     CONTINUE 
          QV(0)=EMBM + EPBM -PISR(0,1)-PISR(0,2)
          QV(1)=-PISR(1,1)+PISR(1,2)
          QV(2)=-PISR(2,1)+PISR(2,2)
          QV(3)=EMBM - EPBM - PISR(3,1)-PISR(3,2)
#if 0
          RSSQ=QV(0)**2-QV(1)**2-QV(2)**2-QV(3)**2
#else
          PA=SQRT(QV(1)*QV(1)+QV(2)*QV(2)+QV(3)*QV(3))
          RSSQ=(QV(0)-PA)*(QV(0)+PA)
#endif
          IF( RSSQ .LE. 0.0 ) THEN
            RS=-SQRT(-RSSQ)
            RETURN
          ENDIF
          RS=SQRT(RSSQ)           
       RETURN
       END

