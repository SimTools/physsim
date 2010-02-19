C***********************************************************************
C*
C*  This programme solves numerically the Schroedinger
C*  equation
C*        dy/dx = z
C*        dz/dx = -mu*(E+iGamma_0/2-V(x))*y
C*  with the running width
C*  and obtain y = y1 and y = y2, whose boundary conditions
C*  are
C*        y1(r0) = r0,   z1(r0) = 1
C*        y2(r0) =  1,   z2(r0) = -beta.
C*  From which Im G(0) is obtained via optical theorem.
C*
C***********************************************************************
 
      SUBROUTINE SGTT0(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                 AMST,VFF,AMSH,BETH,   RSH, SG)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       MODE
      REAL   *8       ALPS, ALP, SN2W, AMSZ, GAMZ, AMSW, AMSB,
     .                AMST, VFF, AMSH, BETH,  RSH, SG
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
      COMPLEX*16      Y1, Y2, Z1, Z2,
     .                DY1, DZ1, DY2, DZ2,
     .                B, BNEW
C--
      PARAMETER( NMAX = 50000 )
#if 1
C>>> For precision plots.
      PARAMETER( EPS = 1.D-11 )
      PARAMETER( TST = 5.D-9  )
#else
C>>> For practical use.
      PARAMETER( EPS = 1.D-6  )
      PARAMETER( TST = 5.D-5  )
#endif
C--
      REAL   *8       XDATA(0:NMAX)
      COMPLEX*16      YDATA(2,0:NMAX),  GL
      COMPLEX*16      GRDATA(0:NMAX)

      DATA AMTS, GAMS / -9999.D0, -9999.D0 /
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical and natural constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xGV2FB = 3.8937966E11
         PI     = ACOS(-1.D0)
         FACT   = 128*PI**2/3 * xGV2FB
         EUL    = 0.5772D0
         IPOT   = 1
         IRGM   = 1
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
C  Initalize effective mass and on-shell width.
C--
10    E    = RSH - 2*AMT
      CALL SETPRM
C--
C  Initalize QCD potential.
C     ALFST = alpha_s(2*m_t)
C--
      IF ( MODE.EQ.0 .OR. ABS(GAM-GAMS)/GAM.GT.0.5D0
     .     .OR. ABS(AMT-AMTS)/AMT.GT.0.2D0 ) THEN
         CALL GTGAMT(ALF,S2W,AMW,AMB,AMT,VTB,GAM)
         CALL SETVPR(ALFS,AMT,AMB,AMZ,GAM,ALFST,ALFSP,RQCD)
         GAMS = GAM
         AMTS = AMT
      END IF
C--
C  Specify boundary conditions for two solutions, y1 and y2,
C  at r = r0.
C--
      R0  = .001D0/AMT
C--
      Y1 = R0
      Z1 = 1
C--
      IF ( IPOT.EQ.1 .OR. IPOT.EQ.2 ) THEN
         Y2 = 1 + BETA*R0*(LOG(LOG(RQCD/R0))
     .                     -(RQCD/R0)*EI(-LOG(RQCD/R0)))
         Z2 = BETA*LOG(LOG(RQCD/R0))
      ELSE IF ( IPOT.EQ.3 ) THEN
         Y2 = 1
         Z2 = -BETA
      ELSE IF ( IPOT.EQ.4 ) THEN
         Y2 = 1
         Z2 = 0
      END IF
C--
C  Consider Higgs effects, if necessary.
C--
      IF ( AMH.LT.1.D3 .AND. BTH.GT.0.D0 ) THEN
         Y2 = Y2 - ZMU*CHIGG*(AMT*BTH)**2*R0*LOG(AMH*R0)
         Z2 = Z2 - ZMU*CHIGG*(AMT*BTH)**2*LOG(AMH*R0)
      END IF
C--
C  Solve Schredinger eq. by Runge-Kutta method.
C--
      NDATA = -1
      B     = 0
      X     = R0
      DXMIN = .01D0/(AMT**3/7000)
#if 0
      XMAX  = 5/SQRT(AMT*GAM)
      DXMAX = XMAX/5000 * MAX(1.D0,SQRT(ABS(E+3)))
#else
      XMAX  = 150/SQRT(AMT*GAM)
      DXMAX = XMAX/150000 * MAX(1.D0,SQRT(ABS(E+3)))
#endif
C--
      DO 2000 N = 0, NMAX
          DX = MIN(MAX(0.002D0*X,DXMIN),DXMAX)
          DO 200 K = 1, 2
             DO 20 I = 1, 10
                CALL  RKUTTA(X,Y1,Z1,DX,DY1,DZ1)
                CALL  RKUTTA(X,Y2,Z2,DX,DY2,DZ2)
                Y1 = Y1 + DY1
                Z1 = Z1 + DZ1
                Y2 = Y2 + DY2
                Z2 = Z2 + DZ2
                X  = X  + DX
20           CONTINUE
C--
C  Save current values of two solutions, y1 and y2.
C--
             NDATA = NDATA + 1
CX           WRITE(20,*) X, Y1, Y2
             XDATA(NDATA)   = X
             YDATA(1,NDATA) = Y1
             YDATA(2,NDATA) = Y2
200       CONTINUE
C--
C  Check if the two solutions match at r = infinity.
C--
          BNEW = - Y2/Y1
          TEST = ABS(B/BNEW-1.)
          IF ( ( TEST.LT.EPS ) .AND.
     .         ( N.GT.500 .OR. X.GT.XMAX ) )     GO TO 250
          B    = BNEW
2000  CONTINUE
C--
C  Iteration over.
C--
      PRINT *,  'Integration steps exceeded',NMAX
C--
C  Two solutions met.
C--
250   CONTINUE
      PRINT *, NDATA, X, B
C--
C  Calculate Im G(0).
C--
      IF ( IRGM.EQ.1 ) THEN
         CAPPA = IMAG(AMT/ZMU)
         QCUT = SQRT(AMT*GAMMA0/2/CAPPA)
         IF ( IPOT.EQ.1 ) THEN
C--
            GI0  = IMAG( (ZMU/4/PI)
     .                   *( 2*QCUT/PI + B + BETA*FX(QCUT*RQCD) ) )
     .             - (ABS(ZMU)**2/4/PI**2)*(GAMMA0/QCUT)
#if 0
                 print *, ' zmu=', zmu, ' qcut=', qcut
                 print *, ' beta=', beta, ' gamma0=', gamma0
                 print *, 'rqcd=', rqcd, ' fx=', fx(qcut*rqcd)
                 print *, 'gi0=', gi0
#endif
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
      ALFT   = 1/ALF
      S      = RSH**2
      DZ     = S/((RSH-AMZ)*(RSH+AMZ))
      SG     = (1/ALFT)**2 / QT**2
     .         *( QT**2 - 2*QT*VE*VT*DZ + (VE**2+AE**2)*VT**2*DZ**2)
     .         *( 1 - 8*ALPS/3/PI )**2
     .         *FACT*GI0/S**2
C--
C  That's it.
C--
      RETURN
      END
