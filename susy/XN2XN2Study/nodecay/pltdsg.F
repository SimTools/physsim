      IMPLICIT    REAL*8  ( A-H, O-Z )
C--
      REAL*8     RS, CS, PHI, SG
      PARAMETER  ( NHEL = 4 )
      INTEGER*4  IHEL(4,NHEL,2)
      DATA IHEL  / -1, +1, -1, -1,
     .             -1, +1, -1, +1,
     .		   -1, +1, +1, -1,
     .             -1, +1, +1, +1,
     .             +1, -1, -1, -1,
     .             +1, -1, -1, +1,
     .		   +1, -1, +1, -1,
     .             +1, -1, +1, +1/
C--
      DATA LOU / 20 /
C--
      REAL*8         RPOL(2), ASFT(3)
C
C========< Entry Point >================================================
C
C--
C  Differential cross section.
C--
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET WINDOW X 2 10 Y 2 8'')')
      WRITE(LOU,'(''SET LIMITS X -1. 1. Y 0. 3.'')')
      WRITE(LOU,'(''(SET LIMITS X -1. 1. Y 0.001 10.'')')
      WRITE(LOU,'(''(SET SCALE Y LOGARITHMIC'')')
      WRITE(LOU,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOU,'(''SET TITLE  SIZE 4.0'')')
      WRITE(LOU,'(''TITLE 5.7 1.2 ''''cos Q'')')
      WRITE(LOU,'(''CASE          ''''    G'')')
#if 0
      WRITE(LOU,'(''TITLE 0.2 3.7 ANGLE 90.'',
     .                            '' ''''dS/dW (fb)'')')
      WRITE(LOU,'(''CASE          '''' G  F     '')')
#else
      WRITE(LOU,'(''TITLE 0.2 3.7 ANGLE 90.'',
     .                            '' ''''dS/dcosQ (fb)'')')
      WRITE(LOU,'(''CASE             '''' G     F     '')')
#endif
C--
C  Standard Model parameters.
C--
      x2PI  = 2*ACOS(-1.D0)
      ALFI  = 128.D0
      ALFS  = 0.12D0
      AMSW  = 80.D0
      AMSZ  = 91.18D0
      AMSH  = 1.D5
      AMST  = 175.D0
      ALF   = 1/ALFI
C--
      IX1   = 2
      IX2   = 2
      AM0   = 500.D0
      AMU   = 400.D0
      AM2   = 250.D0
      TANB  = 3.D0
      AMA   = AM0**2 + AMU**2
      ASFT(1) = 0
      ASFT(2) = 0
      ASFT(3) = 0
C--
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
      CALL INSSCN(ALF,ALFS,AMSW,AMSZ,AM0,AMU,AM2,TANB,AMA,ASFT)
C--
C  Set sqrt(s) and polarization.
C--
      RS   = 500.D0
C      POL  = -0.9D0
      POL  =  0.0D0
C--
      RPOL(1) = (1-POL)/2
      RPOL(2) = (1+POL)/2
      EPSPIN  = 2
C--
C  Angular region for plot.
C     CSMN = cos(theta)_min
C     CSMX = cos(theta)_max
C--
      CSMN = -1
      CSMX = -CSMN
      NCS  = 100
      DCS  = (CSMX-CSMN)/NCS
      PHI  = 0.
C--
C  Loop over cos(theta).
C--
      SGSUM = 0
      WRITE(LOU,'(''SET ORDER X Y '')')
      DO 1000   ICS = 0, NCS
         CS = CSMN + DCS*ICS
         SGSM = 0
         DO 100 IHE = 1, 2
            DO 10 ICASE = 1, 4
               CALL DSGNXX(IX1,IX2,RS,IHEL(1,ICASE,IHE),CS,PHI,SG)
               SGSM = SGSM + RPOL(IHE)*SG/EPSPIN
10          CONTINUE
100      CONTINUE
         SGSUM = SGSUM + SGSM*x2PI*DCS
C--
#if 0
         WRITE(LOU,'(2E15.5)') CS, SGSM
#else
         WRITE(LOU,'(2E15.5)') CS, SGSM*x2PI
#endif
1000  CONTINUE
      WRITE(LOU,'(''JOIN'')')
      WRITE(LOU,'(''( sqrt<s> ='',F8.1,
     .                     '' sigma = '',E15.5,'' fb'')')
     .                     RS, SGSUM
C--
      STOP
      END
