C*
C*  This program plots 1/alpha_i(Q^2) for given alpha, sin^2W,
C*  and alpha_s.
C*
      REAL     *4  B(3,2), AIZ(3), AIQ(3), DAI(3)
      DATA  LOU    / 20 /
      CHARACTER*4  TITLE(2)
      DATA  TITLE  / 'MSSM', 'SM  ' /
C--
C  Constants at Q**2 = MZ**2
C--
      x2PI  = 2*ACOS(-1.)
      x4PI  = 2*x2PI
      ALF   = 1/128.9
      S2W   = 0.23165
      ALFS  = 0.120
      AMZ   = 91.186
C--
      AIZ(1) = (3/5.)*(1-S2W)*(1/ALF)
      AIZ(2) =           S2W *(1/ALF)
      AIZ(3) =                (1/ALFS)
      DAI(1) = 0.11
      DAI(2) = 0.11
      DAI(3) = 0.278
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
      B(1,1) = (3/5.)*( (10/3.)*NG + 1 )
      B(2,1) = 2*NG - 6 + 1
      B(3,1) = 2*NG - 9
C--SM
      B(1,2) = (3/5.)*( (20/9.)*NG + 1/6. )
      B(2,2) = (4/3.)*NG - 22/3. + 1/6.
      B(3,2) = (4/3.)*NG - 11
C--
C  Prepare TOPDRAW data.
C--
      QMN = LOG(1.E2)
      QMX = LOG(1.E17)
      NQ  = 2
      DQ  = (QMX-QMN)/NQ
C--
      DO 1000 I = 1, 2
C--
         WRITE(LOU,'(''(****** '',1A4)') TITLE(I)
         WRITE(LOU,'(''NEW FRAME      '')')
         WRITE(LOU,'(''SET FONT DUPLEX'')')
         WRITE(LOU,'(''SET WINDOW X 2.5 11.5 Y 2 9'')')
C         WRITE(LOU,'(''SET LIMITS X 2 17 Y 0 60'')')
         WRITE(LOU,'(''SET LIMITS X 1.E2 1.E17 Y 0 60'')')
         WRITE(LOU,'(''SET SCALE X LOG'')')
         WRITE(LOU,'(''TITLE 10 8.5 SIZE 4 '''''',1A4)') TITLE(I)
         WRITE(LOU,'(''SET LABELS SIZE 3.0'')')
         WRITE(LOU,'(''TITLE 1.2 5.2 SIZE 4.0 ANGLE 90 ''''1/A'')')
         WRITE(LOU,'(''CASE                            ''''  G'')')
C         WRITE(LOU,'(''TITLE 6.0 1 SIZE 3.5 ''''LOG0101(M0X1) (GeV)'')')
         WRITE(LOU,'(''TITLE 6.0 1 SIZE 3.5 ''''M (GeV)'')')
         WRITE(LOU,'(''CASE                 ''''G      '')')
         WRITE(LOU,'(''SET ORDER X Y '')')
C--
         DO 100 J = 1, 3
            DO 10 IE = -1, 1, 1
               DO 1 IQ = 0, NQ
                  Q  = QMN + IQ*DQ
                  Q  = EXP(Q)
                  AL = (2/x4PI)*LOG(Q/AMZ)
                  AIQ(J) = AIZ(J) + IE*DAI(J) - B(J,I)*AL
C                  WRITE(LOU,'(E12.3,F12.3)') LOG10(Q), AIQ(J)
                  WRITE(LOU,'(E12.3,F12.3)') Q, AIQ(J)
1              CONTINUE
               WRITE(LOU,'(''JOIN'')')
10             CONTINUE
100      CONTINUE
1000  CONTINUE
      PRINT *, '(END)'
C--
C  That's it.
C--
      STOP
      END
