 
      SUBROUTINE  GRQSUM(E,SIGTOT,DLTAFB)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      REAL*8          E, SIGTOT, DLTAFB
C--
C  Parameters.
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
C  Tabulated Green's functions.
C--
      COMMON /STHIRD/ PS, TILG
      REAL   *8       PS(0:100)
      COMPLEX*16      TILG(0:100)
C--
      COMMON /PTHIRD/ PP, TILF
      REAL   *8       PP(0:100)
      COMPLEX*16      TILF(0:100)
C--
C  Coefficients for flat and cos(theta) terms.
C--
      COMMON /COEFF/  TEVEN, TODD1, TODD2, TVTXC, TTOT
      REAL   *8       TEVEN, TODD1, TODD2, TVTXC, TTOT
C--
C  Variables related to optical theorem.
C--
      COMMON /OPTBCD/ B0, B1, BETA, OPT
      REAL   *8       B0, B1
      COMPLEX*16      BETA, OPT
C--
      COMMON /IMG/    B
      COMPLEX*16      B
C--
C  Other variables.
C--
      REAL   *8       T0VC(200), T0FI(200), T1FI(200), T1SP(200),
     .                INTEG0, INTEG1, I0, I1, NORM
C     REAL   *8       INT0WA, I0WA, LNT0WA, LNT2WA
C--
      REAL   *8       LNT0FI, LNT1FI
C
C===========< Entry Point >=============================================
C
C--
C  Calculate cross section using optical theorem.
C--
      TWIMG = 2*DIMAG( EFFMT/4/PI*(OPT+B) )
     .          - ABS(EFFMT)**2/2/PI**2 * GAMMA0(E)/ALAMB
C--
C  Calculate S-wave contributions.
C--
      DO 100 N = 1, 100
         P       = PS(N)
         T0VC(N) = TEVEN * TVTXC * ABS(TILG(N))**2
         T0FI(N) = TEVEN * CF*ALFSP/P * LNT0FI(N)
100   CONTINUE
C--
C  Calculate SP-interference contributions.
C--
      DO 200 N = 1, 100
         P       = PP(N)
         T1SP(N) = TODD1 * P/MT * 2*REAL( TILG(N)*CONJG(TILF(N)) )
200   CONTINUE
C--
C  Calculate contributions from final state interactions.
C--
      DO 300  N = 3, 98
         P       = PS(N)
         T1FI(N) = TODD2*KAPPA * CF*ALFSP/PI/P * LNT1FI(N)
300   CONTINUE
C--
C  Initialize integration step.
C--
      DP    = ALAMB/50
C--
C  Integrate symmetric part.
C--
      INTEG0 = 0
      SIMP1  = 0
      DO 400  N = 1, 99, 2
         SIMP2 =  PS(N)**2*(T0VC(N)+T0FI(N))*RUNWID(E,PS(N))
         SIMP3 =  PS(N+1)**2*(T0VC(N+1)+T0FI(N+1))
     .                      *RUNWID(E,PS(N+1))
         INTEG0 = INTEG0 + ( SIMP1 + 4*SIMP2 + SIMP3 )
         SIMP1 = SIMP3
400   CONTINUE
C--
C  Integrate cos(theta) part.
C--
      INTEG1 = 0
      SIMP1  = 0
      DO 500 N = 1, 97, 2
         SIMP2 = PS(N)**2*(T1SP(N)+T1FI(N))*RUNWID(E,PS(N))
         SIMP3 = PS(N+1)**2*(T1SP(N+1)+T1FI(N+1))
     .                     *RUNWID(E,PS(N+1))
         INTEG1 = INTEG1 + ( SIMP1 + 4*SIMP2 + SIMP3 )
         SIMP1 = SIMP3
500   CONTINUE
C--
C  Loop completed. Now calculate total cross section
C  and FB-asymmetry.
C--
      I0     = INTEG0*DP/6
      I1     = INTEG1*DP/6
      NORM   = 3*ALF**2/2/MT**4
      GV2PB  = 3.8937966D8
      SIGTOT = NORM * I0 * GV2PB
      DLTAFB = I1/2/I0
C>>>
      WRITE(*,*) 'SIGTOT, DELFB = ', REAL(SIGTOT), REAL(DLTAFB)
C>>>
C--
C  That's it.
C--
      RETURN
      END
