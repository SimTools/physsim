CINCLUDE DSGNXX
CINCLUDE SGNXXA
CINCLUDE XNMASS
CINCLUDE SFMASS
CINCLUDE USORTD
C*
C*  Differential cross section for e+e- --> X+ X-.
C*    RS = sqrt(S)
C*
C*  5/16/91  K.Fujii       This version can handle X+1 X-2 case.
C*
      IMPLICIT    REAL*8  ( A-H, O-Y )
      IMPLICIT    COMPLEX*16 ( Z )
C--
      REAL   *8   GZL(2), AMXN(4), AMX(2), GMX(2), AMSF2(7),
     .            AMSE(2), GMSE(2)
      COMPLEX*16  GZX(2), GSEX(2,2)
      REAL   *8   SG(0:6)
      REAL   *8   ON(4,4)
      COMPLEX*16  ETA(4)
C--
      CHARACTER*7 JOIN(1:4)
      DATA JOIN   / 'SOLID  ',
     .              'DOTDASH',
     .              'DASH   ',
     .              'DOT    '/
C--
C     SUSY(1,*,*) = mu
C         (2,*,*) = m_2
C         (3,*,*) = tan(beta)
C         (4,*,*) = m_0
C--
      REAL   *8   SUSY(4,3,9)
      DATA SUSY /
C    .  400.,180.,+50.,200.,  180.,180.,+50.,200.,  180.,400.,+50.,200.,
     . 1000.,400.,-10.,268.,  180.,180.,+50.,200.,  180.,400.,+50.,200.,
C     . 1000.,200.,+10.,500.,  180.,180.,+50.,200.,  180.,400.,+50.,200.,
C    .  400.,250.,+3.,500.,  180.,180.,+50.,200.,  180.,400.,+50.,200.,
C     .  500.,250.,+2.,  70.,  250.,250.,+2.,  70.,  250.,500.,+2.,  70.,
     .  400.,250.,+3., 500.,  250.,250.,+3., 500.,  250.,400.,+3., 500.,
     .  500.,200.,+2., 500.,  200.,200.,+2., 500.,  200.,500.,+2., 500.,
     .  500.,200.,-8., 500.,  200.,200.,-8., 500.,  200.,500.,-8., 500.,
     .  500.,200.,+8., 500.,  200.,200.,+8., 500.,  200.,500.,+8., 500.,
     .  500.,200.,-2.,1000.,  200.,200.,-2.,1000.,  200.,500.,-2.,1000.,
     .  500.,200.,+2.,1000.,  200.,200.,+2.,1000.,  200.,500.,+2.,1000.,
     .  500.,200.,-8.,1000.,  200.,200.,-8.,1000.,  200.,500.,-8.,1000.,
     .  500.,200.,+8.,1000.,  200.,200.,+8.,1000.,  200.,500.,+8.,1000./
      DATA LOU / 20 /
C
C========< Entry Point >================================================
C
C--
C  Differential cross section.
C--
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET WINDOW X 2 10 Y 2 8'')')
      WRITE(LOU,'(''SET LIMITS X -1. 1. Y 0. 0.010'')')
      WRITE(LOU,'(''(SET LIMITS X -1. 1. Y 0.001 10.'')')
      WRITE(LOU,'(''(SET SCALE Y LOGARITHMIC'')')
      WRITE(LOU,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOU,'(''SET TITLE  SIZE 4.0'')')
      WRITE(LOU,'(''TITLE   6.5 7.5 SIZE 2.5'',
     .            '' '''' e2+3e2-3 --> C013S24203C013S24203'')')
      WRITE(LOU,'(''CASE                   '',
     .            '' ''''  X X X X     GUVVMVVX XGUVVMVVX X'')')
      WRITE(LOU,'(''TITLE 5.7 1.2 ''''cos Q'')')
      WRITE(LOU,'(''CASE          ''''    G'')')
      WRITE(LOU,'(''TITLE 0.2 3.7 ANGLE 90.'',
     .                            '' ''''dS/dW (pb)'')')
      WRITE(LOU,'(''CASE          '''' G  F     '')')
      WRITE(LOU,'(''SET ORDER X Y'')')
C--
C  Define constants.
C--
      xPI       = ACOS(-1.)
      x2PI      = 2*xPI
      x4PI      = 2*x2PI
      S2W       = 0.23
      C2W       = 1 - S2W
      SNW       = SQRT(S2W)
      CSW       = SQRT(C2W)
      ALF       = 1.D0/128D0
      SQ2       = SQRT(2.D0)
      AMZ       = 91.17
      AMW       = 80.0
      GMZ       = 2.5
C--
C  Angular region for plot.
C--
      CSMN = -1
      CSMX = -CSMN
      NCS  = 100
      DCS  = (CSMX-CSMN)/NCS
      PHI  = 0.
C--
C  Select reaction type and parameter set.
C--
      IXM       = 1
C      IXM       = 1
      IXP       = 1
C      ISET      = 2
      ISET      = 1
      if (IXM .EQ. IXP) THEN
         FACT = 0.5D0
      ELSE
         FACT = 1.D0
      ENDIF
C--
C  Set sqrt(s) and polarization.
C--
      RS        = 500.
      POL       = 0.
C--
C  Couplings of electron to Z and gamma.
C--
      C0        = SQRT(x4PI*ALF)
      CW        = C0/SNW
      CZ        = CW/CSW
      QE        = -1
      T3LE      = -0.5D0
      T3RE      =  0
C--
      GZL(1)    = CZ*(T3LE-QE*S2W)
      GZL(2)    = CZ*(T3RE-QE*S2W)
C--
      WRITE(LOU,'(''TITLE 7.2 7.0 SIZE 2.5'',
     .            '' ''''2s0O = '',F4.1,'' TeV'')')   RS/1000
      WRITE(LOU,'(''CASE          ''''M UD      '')')
      WRITE(LOU,'(''TITLE 2.5 7.5 SIZE 2.'',
     .            '' ''''m0N013S241   = '',F5.0,'' GeV'')')
     .            SUSY(4,1,ISET)
      WRITE(LOU,'(''CASE '''' XGUVVMVVX   '')')
      WRITE(LOU,'(''TITLE 2.5 7.2 SIZE 2.'',
     .                '' ''''tanB  = '',F5.0)') SUSY(3,1,ISET)
      WRITE(LOU,'(''CASE ''''   G '')')
C--
      WRITE(LOU,'(''TITLE 2.5 6.9 SIZE 2.'',
     .            '' ''''(M,M021) = ('',F5.0,'','',F5.0,''):'')')
     .            SUSY(1,1,ISET), SUSY(2,1,ISET)
      WRITE(LOU,'(''CASE '''' G  X X   '')')
      WRITE(LOU,'(''SET ORDER X Y'')')
      WRITE(LOU,'(''5.7 6.9;6.5,6.9;JOIN TEXT SOLID'')')
      WRITE(LOU,'(''TITLE 2.5 6.6 SIZE 2.'',
     .            '' ''''       = ('',F5.0,'','',F5.0,''):'')')
     .            SUSY(1,2,ISET), SUSY(2,2,ISET)
      WRITE(LOU,'(''SET ORDER X Y'')')
      WRITE(LOU,'(''5.7 6.6;6.5,6.6;JOIN TEXT DOTDASH'')')
      WRITE(LOU,'(''TITLE 2.5 6.3 SIZE 2.'',
     .            '' ''''       = ('',F5.0,'','',F5.0,''):'')')
     .            SUSY(1,3,ISET), SUSY(2,3,ISET)
      WRITE(LOU,'(''SET ORDER X Y'')')
      WRITE(LOU,'(''5.7 6.3;6.5,6.3;JOIN TEXT DASH'')')
C--
C  Loop over SUSY parameters.
C--
      DO 100 ICASE = 1, 3
C--
C  Set SUSY parameters.
C--
         AMU       = SUSY(1,ICASE,ISET)
         AM2       = SUSY(2,ICASE,ISET)
         TNB       = SUSY(3,ICASE,ISET)
         AM0       = SUSY(4,ICASE,ISET)
         BETA      = MOD(ATAN(TNB)+xPI,xPI)
C--
C  Widths and masses.
C--
         CALL XNMASS(AM2,AMU,BETA,S2W,AMZ,AMXN,ON,ETA)
         CALL SFMASS(AM0,AM2,AMU,BETA,ALF,S2W,AMZ,AMSF2)
         AMSE(1) = SQRT(AMSF2(5))
         AMSE(2) = SQRT(AMSF2(7))
C--
         AMX(1)    = AMXN(IXM)
         AMX(2)    = AMXN(IXP)
         GMX(1)    = 0.
         GMX(2)    = 0.
C--
C  Coupling of neutralino to Z.
C--
         ZV        =  CZ/2*(0.D0,1.D0)*IMAG(ETA(IXM)*CONJG(ETA(IXP)))
     .                *( ON(IXM,3)*ON(IXP,3)-ON(IXM,4)*ON(IXP,4) )
         ZA        = -CZ/2            *DBLE(ETA(IXM)*CONJG(ETA(IXP)))
     .                *( ON(IXM,3)*ON(IXP,3)-ON(IXM,4)*ON(IXP,4) )
         IF ( IXM.EQ.IXP ) ZA = ZA/2
C--
         GZX(1)    = ZV - ZA
         GZX(2)    = ZV + ZA
         GSEX(1,1) = SQ2*CZ*ETA(IXM)
     .                *(T3LE*CSW*ON(IXM,2)+(QE-T3LE)*SNW*ON(IXM,1))
         GSEX(2,1) = SQ2*CZ*ETA(IXP)
     .                *(T3LE*CSW*ON(IXP,2)+(QE-T3LE)*SNW*ON(IXP,1))
         GSEX(1,2) = -SQ2*CZ*CONJG(ETA(IXM))*QE*SNW*ON(IXM,1)
         GSEX(2,2) = -SQ2*CZ*CONJG(ETA(IXP))*QE*SNW*ON(IXP,1)
C--
C  Convert GSEX convention to match SUSY manual.
C--
         GSEX(1,1) = CONJG(GSEX(1,1))
         GSEX(2,1) = CONJG(GSEX(2,1))
         GSEX(1,2) = CONJG(GSEX(1,2))
         GSEX(2,2) = CONJG(GSEX(2,2))
C--
         PRINT *, ' GZX       = ', GZX
         PRINT *, ' GSEX(*,1) = ', GSEX(1,1), GSEX(2,1)
         PRINT *, ' GSEX(*,2) = ', GSEX(1,2), GSEX(2,2)
C--
         SGSM = 0
         WRITE(LOU,'(''(                                    '')')
         WRITE(LOU,'(''(  **********************************'')')
         WRITE(LOU,'(''(    AMU   = '',F10.3)') SUSY(1,ICASE,ISET)
         WRITE(LOU,'(''(    AM2   = '',F10.3)') SUSY(2,ICASE,ISET)
         WRITE(LOU,'(''(    TANB  = '',F10.3)') SUSY(3,ICASE,ISET)
         WRITE(LOU,'(''(    AM0   = '',F10.3)') SUSY(4,ICASE,ISET)
         WRITE(LOU,'(''(   ---------------------------------'')')
         WRITE(LOU,'(''(    AM1   = '',F10.3)') AMX(1)
         WRITE(LOU,'(''(    AM2   = '',F10.3)') AMX(2)
         WRITE(LOU,'(''(    AMSEL = '',F10.3)') AMSE(1)
         WRITE(LOU,'(''(    AMSER = '',F10.3)') AMSE(2)
         WRITE(LOU,'(''(  **********************************'')')
         WRITE(LOU,'(''(                                    '')')
C--
         DO 10   ICS = 0, NCS
            CS = CSMN + DCS*ICS
            IF ( ICS.EQ.0 .OR. ICS.EQ.NCS ) THEN
               WGT = 1./3
            ELSE IF ( MOD(ICS,2).EQ.0 ) THEN
               WGT = 2./3
            ELSE
               WGT = 4./3
            ENDIF
            CALL DSGNXX(GZL,GZX,GSEX,AMSE,GMSE,AMZ,GMZ,AMX,GMX,
     .                  RS,POL,CS,PHI,SG)
            WRITE(LOU,'(2E15.5)') CS, SG(0)*FACT
            SGSM = SGSM + SG(0)*WGT*FACT
10       CONTINUE
         WRITE(LOU,'(''JOIN '',A7)') JOIN(ICASE)
         WRITE(LOU,'(''( sqrt<s> ='',F8.1,
     .                  '' sigma = '',E15.5,'' pb'')')
     .                  RS, SGSM*x2PI*DCS
C--
C  Compare numerical integral and analytic integral.
C--
         WRITE(6,'(''( sqrt<s> ='',F8.1,
     .                  '' sig_n = '',E15.5,'' pb'')')
     .                  RS, SGSM*x2PI*DCS
C--
         CALL SGNXXA(GZL,GZX,GSEX,AMSE,GMSE,AMZ,GMZ,AMX,GMX,
     .                  RS,POL,SG)
         WRITE(6,'(''( sqrt<s> ='',F8.1,
     .                  '' sig_a = '',E15.5,'' pb'')')
     .                  RS, SG(0)*FACT
100   CONTINUE
C--
      STOP
      END
