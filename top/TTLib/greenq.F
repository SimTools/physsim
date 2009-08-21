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
C*        y2(r0) =  1,   z2(r0) = -beta.
C*  From which Im G(0) is obtained via optical theorem.
C*
C***********************************************************************
 
      SUBROUTINE GREENQ(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                  AMST,VFF,AMSH,BETH,   RSH,Q, GQ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       MODE
      REAL   *8       ALPS, ALP, SN2W, AMSZ, GAMZ, AMSW, AMSB,
     .                AMST, VFF, AMSH, BETH,  RSH, Q
      COMPLEX*16      GQ
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
#if 0
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
#if 1
      PARAMETER       ( NQ = 500 )
#else
      PARAMETER       ( NQ = 5000 )
#endif
      REAL   *8       QDATA(0:NQ)
      COMPLEX*16      GRDATA(0:NMAX), GQDATA(0:NQ)

      DATA AMTS, GAMS / -9999.D0, -9999.D0 /
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical and natural constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI    = ACOS(-1.D0)
         IPOT  = 1
         IRGM  = 1
C--
         ALF   = ALP
         S2W   = SN2W
         AMB   = AMSB
         AMW   = AMSW
         AMZ   = AMSZ
         GMZ   = GAMZ
         CHIGG = 1/(4*PI)/(AMW**2*S2W/PI/ALF)
      END IF
C--
C  Check if initialization necessary.
C     MODE = 0 : first call or if alpha_s changed.
C          = 1 : if alpha_s unchanged.
C          = 2 : if only RSH changed.
C          = 3 : if all the parameters unchanged.
C--
      IF ( MODE.EQ.2 )                           GO TO 10
      IF ( MODE.EQ.3 )                           GO TO 550
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
C  Calculate Green's function in configuration space.
C--
      ANRMMX = 0
      DO 300 I = 0, NDATA
          GL     = YDATA(2,I) + B*YDATA(1,I)
          ANRM   = ABS(GL)**2
          ANRMMX = MAX( ANRM, ANRMMX )
          GRDATA(I) = GL
          IF ( ANRM .LT. TST*ANRMMX )            GO TO 350
300   CONTINUE
C--
C  Error.
C--
      PRINT *, 'G>(r) does not become sufficiently small.'
C--
C  Calculate Green's function in momentum space.
C--
350   KMAX  = I/2
      QMAX  = SQRT( (3*AMT+E-AMW)*(3*AMT+E+AMW)
     .              *(AMT+E-AMW)*(AMT+E+AMW) )/2/(2*AMT+E)
      DQ    = QMAX/NQ
      DO 500 IQ = 0, NQ
         QDATA(IQ)  = IQ*DQ
         GQDATA(IQ) = 0
         IF ( IQ.EQ.0 ) THEN
            DO 40 K = 0, KMAX-1
               I1 = 2*K
               I2 = I1 + 1
               I3 = I1 + 2
               DR = XDATA(I2) - XDATA(I1)
               GQDATA(IQ) = GQDATA(IQ) +   XDATA(I1)*GRDATA(I1)*DR
     .                                 + 4*XDATA(I2)*GRDATA(I2)*DR
     .                                 +   XDATA(I3)*GRDATA(I3)*DR
40          CONTINUE
            GQDATA(IQ) = - GQDATA(IQ)*ZMU/3
         ELSE
            DO 50 K = 0, KMAX-1
               I1 = 2*K
               I2 = I1 + 1
               I3 = I1 + 2
               DR = XDATA(I2) - XDATA(I1)
               THETA1 = QDATA(IQ)*XDATA(I1)
               THETA2 = QDATA(IQ)*XDATA(I2)
               THETA3 = QDATA(IQ)*XDATA(I3)
               GQDATA(IQ) = GQDATA(IQ) +   SIN(THETA1)*GRDATA(I1)*DR
     .                                 + 4*SIN(THETA2)*GRDATA(I2)*DR
     .                                 +   SIN(THETA3)*GRDATA(I3)*DR
50          CONTINUE
            GQDATA(IQ) = -GQDATA(IQ)*ZMU/QDATA(IQ)/3
         END IF
500   CONTINUE
C--
C  Get Green's function in momentum space from table.
C--
550   IQ = Q/DQ
      IF ( IQ.GE.NQ-1 ) THEN
         GQ  = 0
      ELSE IF ( Q.LE.0.D0 ) THEN
         GQ  = GQDATA(0)
      ELSE
         A1 = Q - QDATA(IQ)
         A2 = QDATA(IQ+1) - Q
         GQ = (A2*GQDATA(IQ)+A1*GQDATA(IQ+1))/DQ
      END IF
#if 0
C     Free Green's function
#if 0
      s1 = amt**2
      e1 = sqrt(s1 + q**2)
      e2 = rsh - e1
      s2 = (e2 - q2)*(e2 + q2)
      gq = 1/((s1 + s2)/(2*amt) - amt + (0.D0, 1.D0)*gam)
#else
      gq = 1/(e - q**2/amt + (0.D0, 1.D0)*gam)
#endif
#endif
C--
C  That's it.
C--
      RETURN
      END
