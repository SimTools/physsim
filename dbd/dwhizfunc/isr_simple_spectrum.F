
      SUBROUTINE ISR_SIMPLE_SPECTRUM(Z1, Z2, xALF0, 
     >       QED, RS, EMBM, EPBM, QV)
C =======================================================================
C
C(Function)
C  Generate ISR spectrum with simple formula
C(Input)
C  xALF0 : alpha0 ( = 1/137.....)
C  Z1    : random number for brems photon
C  Z2    : random number to decide e- or e+ to brems
C  EMBM  : e- beam energh after beamstrahlung
C  EPBM  : e+ beam energy after beamstrahlung
C  RS    : Root-S of initial system
C(Output)
C  EMBM  : e- energy after ISR
C  EPBM  : e+ energy after ISR
C  QV(0:3) : 4 momentum of hard scattering ststem after ISR
C  RS    : Root-S of final system excluding ISR
C  QED   : QED factor
C(Author)
C  A.Miyamoto  21-June-2011 ( extracted from functtf.F )
C =======================================================================
      IMPLICIT NONE
      REAL*8   xALF0, QV(0:3), RS, EMBM, EPBM
      REAL*8   Z1, Z2
#include "SMCONS.inc"
#include "SMPTAB.inc"
      REAL*8   BTE, QED, ZG, AME
      INTEGER*4 IZE

      AME = AMF(1,2,1)

      BTE = (2*xALF0/xPI)*(2*LOG(RS/AME)-1)
      QED = ( 1 + 3*BTE/4 )*(1 + 2*xALF0/xPI*(xPI**2/6.0D0-1/4.D0))
      ZG  = Z1**(1/BTE)
      RS  = SQRT(1-ZG)*RS
      IZE = Z2*1.E6
      IF ( MOD(IZE,2).EQ.0 ) THEN
         EMBM = (1-ZG)*EMBM
      ELSE
         EPBM = (1-ZG)*EPBM
      ENDIF
      QV(0)=EMBM+EPBM
      QV(1)=0
      QV(2)=0
      QV(3)=EMBM-EPBM
      RS   = 2*SQRT(EMBM*EPBM)

      RETURN 
      END

