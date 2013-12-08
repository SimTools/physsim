C***********************************************************************
C*
C*  This programme solves numerically the Schroedinger
C*  equation
C*        dy/dx = z
C*        dz/dx = -mu*(E+iGamma_0/2-V(x))*y
C*  with the RUNNING WIDTH
C*  and obtain y = y1 and y = y2, whose boundary conditions
C*  are
C*        y1(r0) = r0,   z1(r0) = 1
C*        y2(r0) =  1,   z2(r0) = ??? depends on potentials
C*  From which Im G(0) is obtained via optical theorem.
C*
C*  93/06/17  K.Fujii                  This version uses Higgs effects
C*                                     calculated perturbatively
C*                                     instead of those calculated with
C*                                     potential approximation.
C*
C***********************************************************************
 
      SUBROUTINE SGTTHR(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                  AMST,VFF,AMSH,BETH,   RSH,SG)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4 MODE
      REAL*8    ALPS, ALP, SN2W, AMSZ, GAMZ, AMSW, AMSB,
     .          AMST, VFF, AMSH, BETH,  RSH, SG
C--
C  Potential selection.
C     IPOT = 1: QCD( Hagiwara   )
C          = 2: Richardson
C          = 3: Cornell
C          = 4: Martin
C  Width selection.
C     IRGM = 1: running width
C          = 0: constant width
C--
      COMMON /KFFLAG/ IPOT, IRGM
      INTEGER*4       IPOT, IRGM
C--
      COMMON /EFFMAS/ ZMU, BETA
      COMPLEX*16      ZMU, BETA
      COMMON /PARAM1/ E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      REAL   *8       E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
C--
      COMPLEX*16      G1(2), G2(2), GPS(2), B
      DATA AMTS, GAMS / -9999.D0, -9999.D0 /
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical and natural constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xGV2PB = 3.8937966E8
         PI     = ACOS(-1.D0)
         FACT   = 128*PI**2/3 * xGV2PB
         EUL    = 0.5772D0
         IPOT   = 1
         IRGM   = 1
         CALL INITRK
C--
         ALF   = ALP
         S2W   = SN2W
         AMB   = AMSB
         AMW   = AMSW
         AMZ   = AMSZ
         GMZ   = GAMZ
         CHIGG = 1/(4*PI)/(AMW**2*S2W/PI/ALF)
C--
         SCW   = SQRT(S2W*(1-S2W))
         QE    = -1
         VE    = ( -1.D0/2 - 2*QE*S2W )/2/SCW
         AE    = ( -1.D0/2            )/2/SCW
         QT    = 2.D0/3
         VT    = ( +1.D0/2 - 2*QT*S2W )/2/SCW
      END IF
C--
C  Check if initialization necessary.
C     MODE = 0 : first call or if alpha_s changed.
C          = 1 : if alpha_s unchanged.
C          = 2 : if only RSH changed.
C--
      IF ( MODE.EQ.2 )                           GO TO 10
C--
C  Initialize variable prameters in /PARAM1/.
C--
      ALFS = ALPS
      AMT  = AMST
      VTB  = VFF
      AMH  = AMSH
      BTH  = BETH
C--
C  Calculate energy measured from threshold.
C--
10    E    = RSH - 2*AMT
C--
C  Initalize QCD potential.
C     ALFST = alpha_s(2*m_t)
C--
      CALL GTGAMT(ALF,S2W,AMW,AMB,AMT,VTB,GAM)
      IF ( MODE.EQ.0 .OR. ABS(GAM-GAMS)/GAM.GT.0.5D0
     .     .OR. ABS(AMT-AMTS)/AMT.GT.0.2D0 ) THEN
         CALL SETVPR(ALFS,AMT,AMB,AMZ,GAM,ALFST,ALFSP,RQCD)
         GAMS = GAM
         AMTS = AMT
         PRINT *, ' m_t = ', AMT, ' Gamma_t w/o QCD = ', GAM
         PRINT *, ' alpha_s(mz     ) = ', ALFS
         PRINT *, ' alpha_s(mt     ) = ', ALFST
         PRINT *, ' alpha_s(mt*alfs) = ', ALFSP
      END IF
C--
C  Initalize effective mass and on-shell width.
C--
      CALL SETPRM
C--
C  Specify boundary conditions for two solutions, y1 and y2,
C  at r = r0.
C--
      R0  = .001/AMT
C--
      G1(1) = R0
      G1(2) = 1
      IF ( IPOT.EQ.1 .OR. IPOT.EQ.2 ) THEN
         G2(1) = 1 + BETA*R0*(LOG(LOG(RQCD/R0))
     .                     -(RQCD/R0)*EI(-LOG(RQCD/R0)))
         G2(2) = BETA*LOG(LOG(RQCD/R0))
      ELSE IF ( IPOT.EQ.3 ) THEN
         G2(1) = 1
         G2(2) = -BETA
      ELSE IF ( IPOT.EQ.4 ) THEN
         G2(1) = 1
         G2(2) = 0
      END IF
C--
C  Consider Higgs effects, if necessary.
C--
      IF ( AMH.LT.1.D3 .AND. BTH.GT.0.D0 ) THEN
         CALL FHIGGS(ALF,S2W,AMW,AMT,AMH,BTH,FACH)
      ELSE
         FACH = 1
      END IF
C--
C  Solve Schredinger eq. by Runge-Kutta method.
C--
#if defined(NARROW_WIDTH)
C>>>
C For narrow width resonances.
      REND  = 15/SQRT(AMT*GAM) * MAX(1.D0,(2+E)/1.5D0)
C>>> 20090815 KF
#else
C     REND  = 20/SQRT(AMT*GAM)
      REND  = 150/SQRT(AMT*GAM)
C<<< 20090815
#endif
C>>>
      CALL RHS(R0,G1,GPS)
      CALL ODEINT(G1,GPS,R0,REND)
      CALL RHS(R0,G2,GPS)
      CALL ODEINT(G2,GPS,R0,REND)
      B = - G2(1)/G1(1)
C--
C  Calculate Im G(0).
C--
      IF ( IRGM.EQ.1 ) THEN
         CAPPA = IMAG(AMT/ZMU)
         QCUT = SQRT(AMT*GAMMA0/2/CAPPA)
CX       QCUT = SQRT( (3*AMT+E-AMW)*(3*AMT+E+AMW)
CX   .               *(AMT+E-AMW)*(AMT+E+AMW) )/2/(2*AMT+E)
         IF ( IPOT.EQ.1 ) THEN
C--
            GI0  = IMAG( (ZMU/4/PI)
     .                   *( 2*QCUT/PI + B + BETA*FX(QCUT*RQCD) ) )
     .             - (ABS(ZMU)**2/4/PI**2)*(GAMMA0/QCUT)
         ELSE IF ( IPOT.EQ.2 ) THEN
            GI0  = IMAG( (ZMU/4/PI)
     .                   *( 2*QCUT/PI + B + BETA*F1(QCUT*RQCD) ) )
     .             - (ABS(ZMU)**2/4/PI**2)*(GAMMA0/QCUT)
         ELSE IF ( IPOT.EQ.3 ) THEN
            GI0  = IMAG( (ZMU/4/PI)*( 2*QCUT/PI + B
     .                     + BETA*( EUL - 1 + LOG(QCUT*R0) ) ) )
     .             - (ABS(ZMU)**2/4/PI**2)*(GAMMA0/QCUT)
         ELSE IF ( IPOT.EQ.4 ) THEN
            GI0  = IMAG( (ZMU/4/PI)*( 2*QCUT/PI + B ) )
     .             - (ABS(ZMU)**2/4/PI**2)*(GAMMA0/QCUT)
         END IF
      ELSE
         GI0  = IMAG( (ZMU/4/PI)*B )
      END IF
C--
C  Calculate cross section.
C--
CO       S      = RSH**2
CO       SG     = (1/ALFT)**2 *( 1 + 2.25D0*((1-2.6667D0*S2W)
CO      .                                    /(16*S2W*(1-S2W)))**2
CO      .                                  *(S/(RSH-AMZ)/(RSH+AMZ))**2 )
CO      .         *( 1 - 8*ALFST/3/PI )**2
CO      .         *FACT*GI0/S**2
C--
CX    ALFT   = 1/ALF - 1.415*DLOG(RSH/AMZ)
      ALFT   = 1/ALF
      S      = RSH**2
      DZ     = S/((RSH-AMZ)*(RSH+AMZ))
      SG     = FACH*(1/ALFT)**2 / QT**2
     .         *( QT**2 - 2*QT*VE*VT*DZ + (VE**2+AE**2)*VT**2*DZ**2)
     .         *( 1 - 8*ALFST/3/PI )**2
     .         *FACT*GI0/S**2
C--
C  That's it.
C--
350   RETURN
      END
C*
C* (Purpose)
C*     This program computes the right hand side of the differential
C*     equation:    G'(I) = GPR(I).
 
      SUBROUTINE RHS(R,G,GPR)
 
      IMPLICIT REAL*8 ( A-F, P-Z )
C--
      COMMON /PARAM1/ E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      REAL   *8       E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      COMMON /EFFMAS/ ZMU, BETA
      COMPLEX*16      ZMU, BETA
C--
      REAL   *8       R
      COMPLEX*16      G(2), GPR(2)
      COMPLEX*16      ZI
      PARAMETER     ( ZI = (0.D0,1.D0) )
C
C========< Entry Point >================================================
C
C--
C  Calculate right-hand sides.
C--
      GPR(1) = G(2)
      GPR(2) = ZMU*(V(R)-E-ZI*GAMMA0/2)*G(1)
C--
C  That's it.
C--
      RETURN
      END
