CC**********************************************************************
C*
C*  SUBROUTINE GTCEFF(E)
C*
C*  This subroutine calculates the coefficients, t_even, t_odd1, t_odd2,
C*  and t_vtxc necessary for the calculation of T_0, T_1 and T_2.
C*
C*       t_even = ( -2/3 + chi*geV*gtV )**2 + (chi*geA*gtV)**2
C*
C*       t_odd1 =  -2/3*chi*geA*gtA + 2*chi**2*geV*geA*gtV*gtA
C*
C*       t_odd2 = (-2*chi*geA*gtV)*( -2/3 + chi*geV*gtV )
C*
C*       t_vtxc = 1 - 16*alpha_s/3/pi
C*                  - CF*alpha_s/2/pi* 2*h_bWg(YCUT,r)
C*
C*  Using these coefficients, T_0 and T_1 are given by
C*
C*       T0 = T0_VC + T0_FI
C*       T1 = T1_FI + T1_SP
C*
C*  Additional contribution comes from the 'wrong assignment' diagram:
C*
C*       T0_WA and T2_WA
C*
C*  Here,
C*
C*       T0_VC = t_even*t_vtxc * !tilde{G}(p;E)!^2
C*       T0_FI = t_even * CF*4*pi*alpha_s *
C*         \int d^3{q}/(2pi)^3  1/!p_g!^3
C*               * 2*Im[ conj(tilde{G}(q;E))*tilde{G}(p;E) ] * pi/2
C*
C*       T1_SP = t_odd1*p/m_t*2*Re[ tilde{G}(p;E)*conj(tilde{F}(p;E)) ]
C*       T1_FI = t_odd2*kappa*CF*4*pi*alpha_s *
C*         \int d^3{q}/(2pi)^3  1/!p_g!^3
C*               2*Re[ conj(tilde{G}(q;E))*tilde{G}(p;E) ] * P1
C*
C*       T0_WA = t_even * CF*4*pi*alpha_s *
C*         \int d^3{q}/(2pi)^3  1/!p_g!^3 * !tilde{G}(q;E)!^2
C*           *( Phi00 - kappa^2*Phi55 + 2/3*P0*kappa^2*Phi55 )
C*       T2_WA = t_even * kappa^2 * CF*4*pi*alpha_s *
C*         \int d^3{q}/(2pi)^3  1/!p_g!^3 * !tilde{G}(q;E)!^2
C*           *(-2/3*P2)*Phi55
C*
C*  Inputs:
C*       alft : alpha_s(m_t) ; stored in  common/QCDpar/...
C*       mt, mZ, pi;  stored in  common/TTPARM/.
C*  Outputs:
C*       t_even, t_odd1, t_odd2, t_vtxc; stored in  common/coeff/...
C*
C* (Update Record)
C*    93/05/08  K.Fujii                New argument E is added to
C*                                     take account of E-dependence
C*                                     of the Z/gamma propagator
C*                                     factor.
C*
CC**********************************************************************
 
      SUBROUTINE GTCEFF(E)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
C--
      REAL   *8       E
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /QCDPAR/ ALFST, ALFSP, RQCD
      REAL*8          ALFST, ALFSP, RQCD
C--
      COMMON /COEFF/  TEVEN, TODD1, TODD2, TVTXC, TTOT
      REAL   *8       TEVEN, TODD1, TODD2, TVTXC, TTOT
C--
      REAL   *8       CHI, GEV, GEA, GTV, GTA
C--
      REAL   *8       HBWG
C
C========< Entry Point >================================================
C
C--
C  Initialize parameters.
C--
      S   = (2*MT+E)**2
      CHI = 1/( 4*S2W*(1-S2W) ) * S/( S - MZ**2 )
      FAC = ( 4*MT**2/S )**2
C--
      GEV = -0.5D0 + 2.D0*S2W
      GEA = +0.5D0
      GTV = +0.5D0 - 4.D0/3.D0*S2W
      GTA = -0.5D0
C--
      TEVEN = (( -2.D0/3.D0 + CHI*GEV*GTV )**2 + (CHI*GEA*GTV)**2)*FAC
      TODD1 = (-2.D0/3.D0*CHI*GEA*GTA + 2*CHI**2*GEV*GEA*GTV*GTA)*FAC
      TODD2 = (-2*CHI*GEA*GTV)*( - 2.D0/3.D0 + CHI*GEV*GTV )*FAC
C--
      TVTXC = 1 - 16*ALFST/3/PI
     .            - CF*ALFST/2/PI* 2*HBWG(RWT,YCUT)
C--
C  That's it.
C--
      RETURN
      END
