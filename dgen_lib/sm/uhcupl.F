C***********************************************************************
C*
C*==============================----------------------------------==
C* Subroutine UHCUPL(ALPHA,SW2, GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1,
C*                   GW,GWWZ,GWWA)
C*===================------------===================================
C*
C* (Ppurpose)
C*    This subroutine sets up the coupling constants in the SM.
C*    The array of couplings specifies the chirality of the INCOMING
C*    fermion. G?(1) denotes a left-handed coupling, and G?(2) a right-
C*    handed coupling.
C* (Inputs)
C*      ALPHA    : (R*8) : fine structure constant.
C*      SW2      : (R*8) : square of sine of Weinberg angle
C* (Outputs)
C*      GAL(2)   : (R*8) : coupling with A of charged leptons.
C*      GAU(2)   : (R*8) : coupling with A of up-type quarks.
C*      GAD(2)   : (R*8) : coupling with A of down-type quarks.
C*      GWF(2)   : (R*8) : coupling with W-,W+ of fermions.
C*      GZN(2)   : (R*8) : coupling with Z of neutrinos.
C*      GZL(2)   : (R*8) : coupling with Z of charged leptons.
C*      GZU(2)   : (R*8) : coupling with Z of up-type quarks.
C*      GZD(2)   : (R*8) : coupling with Z of down-type quarks.
C*      G1(2)    : (R*8) : unit coupling of fermions.
C*      GW       : (R*8) : weak coupling constant.
C*      GWWZ     : (R*8) : coupling of W-,W+,Z.
C*      GWWA     : (R*8) : coupling of W-,W+,A.
C* (Relation)
C*    Calls no subroutines.
C* (Update Record)
C*    6/22/90  K.Fujii       Derived from COUPLx in HELAS.
C*    6/14/92  K.Fujii       New HELAS convension which differs from
C*                           the old in the signs of ffV vertices.
C*
C***********************************************************************
 
      SUBROUTINE UHCUPL(ALPHA,SW2, GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1,
     .                  GW,GWWZ,GWWA)
 
      REAL*8    ALPHA,SW2,GW,GWWZ,GWWA
      REAL*8    GAL(2),GAU(2),GAD(2),GWF(2),GZN(2),GZL(2),
     .          GZU(2),GZD(2),G1(2)
C
C========< Entry Point >================================================
C
C--
C  Define some constants.
C--
      FOURPI = 4*ACOS(-1.D0)
      E      = SQRT(ALPHA *FOURPI)
      SW     = SQRT(SW2)
      CW     = SQRT(1-SW2)
      EZ     = E/(SW*CW)
C--
C  Calculate coupling constants involving fermions.
C--
      GAL(1) = E
      GAL(2) = E
      GAU(1) =-E*2.D0/3
      GAU(2) =-E*2.D0/3
      GAD(1) = E/3
      GAD(2) = E/3
      GWF(1) =-E/SQRT(2*SW2)
      GWF(2) = 0
      GZN(1) =-EZ/2
      GZN(2) = 0
      GZL(1) =-EZ*(-1.D0/2+SW2)
      GZL(2) =-EZ*SW2
      GZU(1) =-EZ*(1.D0/2-2.D0/3*SW2)
      GZU(2) = EZ*2.D0/3*SW2
      GZD(1) =-EZ*(-1.D0/2+1.D0/3*SW2)
      GZD(2) =-EZ/3*SW2
      G1(1)  = 1
      G1(2)  = 1
C--
C  Calculate gauge boson self-coupling constants.
C--
      GW     = E/SW
      GWWZ   = E/SW*CW
      GWWA   = E
C--
C  That's it.
C--
      RETURN
      END
