C*
C*  This program calculates alpha_s(M_Z) from alpha and S2W
C*  assuming GUT.
C*
      REAL     *4  B1(2), B2(2), B3(2)
      CHARACTER*4  TITLE(2)
      DATA  TITLE  / 'MSSM', 'SM  ' /
C--
C  Constants at Q**2 = MZ**2
C--
      ALF   = 1/128.8
      x2PI  = 2*ACOS(-1.)
      x4PI  = 2*x2PI
      S2W   = 0.233
      AMZ   = 91.1
C--
C  Beta functions.
C--
      NG    = 3
C--
      PRINT *, ' 1/ALF(MZ) = ', 1/ALF
      PRINT *, ' SIN2W     = ', S2W
      PRINT *, ' AMZ       = ', AMZ, ' GeV'
      PRINT *, ' NG        = ', NG
      PRINT *, ' --- '
C--MSSM
      B1(1) = (10/3.)*NG + 1
      B2(1) = 2*NG - 6 + 1
      B3(1) = 2*NG - 9
C--SM
      B1(2) = (20/9.)*NG + 1/6.
      B2(2) = (4/3.)*NG - 22/3. + 1/6.
      B3(2) = (4/3.)*NG - 11
C--
C  Loop over models.
C--
      DO 10 I = 1, 2
C--
C  MX.
C--
         ALN   = (x4PI/ALF)*(3-8*S2W)/(3*B1(I)-5*B2(I))
         AMX   = AMZ*EXP(ALN/2)
         ALFXI = S2W/ALF - (B2(I)/x4PI)*ALN
         ALFSI = ALFXI   + (B3(I)/x4PI)*ALN
C--
C  Print out results.
C--
         PRINT *, TITLE(I)
         PRINT *, ' MX   = ', AMX, ' GeV'
         PRINT *, ' ALFX = ', 1/ALFXI
         PRINT *, ' ALFS = ', 1/ALFSI
10    CONTINUE
C--
C  That's it.
C--
      STOP
      END
