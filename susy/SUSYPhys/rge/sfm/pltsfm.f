C*
C*  This program plots scalar masses as functions of gluino mass.
C*
C*
      PARAMETER  ( NG = 20, NS = 3 )
      REAL*4  ALFT(3), B(3), F(3), C(2,7), XYDATA(2,0:NG,7,0:NS)
      CHARACTER*7 JOIN(7)
      DATA JOIN / 'SOLID  ', 'SOLID  ', 'DOTDASH', 'DOTDASH',
     .            'DASH   ', 'DASH   ', 'DOT    '/
      DATA LOU  / 20 /
C--
C  Constants at Q**2 = MZ**2
C--
      ALF   = 1/128.8
      x2PI  = 2*ACOS(-1.)
      x4PI  = 2*x2PI
      S2W   = 0.233
      AMZ   = 91.1
      NGN   = 3
C--
C  Beta functions.
C--
C--MSSM
      B(1)  = (10/3.)*NGN + 1
      B(2)  = 2*NGN - 6 + 1
      B(3)  = 2*NGN - 9
C--
C  MX.
C--
      TEE   = (x4PI/ALF)*(3-8*S2W)/(3*B(1)-5*B(2))
      AMX   = AMZ*EXP(TEE/2)
      ALFXI = S2W/ALF - (B(2)/x4PI)*TEE
      ALFSI = ALFXI   + (B(3)/x4PI)*TEE
C--
      PRINT *, ' 1/ALF(MZ) = ', 1/ALF
      PRINT *, ' SIN2W     = ', S2W
      PRINT *, ' AMZ       = ', AMZ, ' GeV'
      PRINT *, ' NGN       = ', NGN
      PRINT *, ' AMX       = ', AMX
      PRINT *, ' --- '
C--
C  Alpha_tilde.
C--
      ALFT(3) = 1/ALFXI/x4PI
      ALFT(2) = ALFT(3)
      ALFT(1) = 3/5.*ALFT(2)
C--
C  f_i.
C--
      DO 10 I = 1, 3
         F(I) = 1/B(I)/ALFT(I) * ( 1 - 1/(1+B(I)*ALFT(I)*TEE)**2 )
10    CONTINUE
C--
C  M_i**2 = M_0**2 + C(1,i)*M**2 + C(2,i)*cos(2*beta)*M_Z**2.
C--
      C(1,1) = 2*( 4/3.*ALFT(3)*F(3) + 3/4.*ALFT(2)*F(2)
     .                               + 1/36.*ALFT(1)*F(1) )
      C(1,2) = C(1,1)
      C(1,3) = 2*( 4/3.*ALFT(3)*F(3) + 4/9.*ALFT(1)*F(1) )
      C(1,4) = 2*( 4/3.*ALFT(3)*F(3) + 1/9.*ALFT(1)*F(1) )
      C(1,5) = 2*( 3/4.*ALFT(2)*F(2) + 1/4.*ALFT(1)*F(1) )
      C(1,6) = C(1,5)
      C(1,7) = 2*ALFT(1)*F(1)
C--
      DO 66 I = 1, 7
66    PRINT *, ' C(1,', I, ') = ', C(1,I)
C--
      C1QSM  = C(1,1) + C(1,2) + C(1,3) + C(1,4)
      C1LSM  = C(1,5) + C(1,6) + C(1,7)
C--
      C(2,1) = -0.5 + 2/3.*S2W
      C(2,2) =  0.5 - 1/3.*S2W
      C(2,3) =      - 2/3.*S2W
      C(2,4) =      + 1/3.*S2W
      C(2,5) =  0.5 -      S2W
      C(2,6) = -0.5
      C(2,7) =      +      S2W
C--
C  Loop over gluino masses.
C--
      AMGMN = 500
      AMGMX = 2000
      DMG   = (AMGMX-AMGMN)/NG
      TNBMN = 1.
      TNBMX = 4.
      DTNB  = (TNBMX-TNBMN)/NS
      DO 2000 IS = 0, NS
         TNB   = TNBMN + DTNB*IS
         AMZCS = AMZ*AMZ*(TNB-1)*(TNB+1)/(TNB*TNB+1)
         DO 200  IG = 0, NG
            AMG  = AMGMN + DMG*IG
            AM   = ALFSI/ALFXI *AMG
            DO 20  I = 1, 7
               XYDATA(1,IG,I,IS) = AMG
               XYDATA(2,IG,I,IS) = AM*AM*C(1,I) + C(2,I)*AMZCS
20          CONTINUE
200      CONTINUE
2000  CONTINUE
C--
C  Print out results.
C--
      DO 3000 IS = 0, NS
         TNB = TNBMN + DTNB*IS
         WRITE(LOU,'(''NEWFRAME'')')
         WRITE(LOU,'(''SET FONT DUPLEX'')')
         WRITE(LOU,'(''( TNB = '',F5.0)') TNB
         WRITE(LOU,'(''SET LIMITS X 500 2000'')')
         WRITE(LOU,'(''SET LIMITS Y 50 3000'')')
         WRITE(LOU,'(''SET SCALE Y LOG'')')
C        WRITE(LOU,'(''SET ORDER X Y'')')
         DO 300 I = 1, 7
            WRITE(LOU,'(''( I = '',I3)') I
            DO 30 IG = 0, NG
               X = XYDATA(1,IG,I,IS)
               Y = XYDATA(2,IG,I,IS)
               WRITE(LOU,*) X, SQRT(Y)
30          CONTINUE
            WRITE(LOU,'(''JOIN '',A7)') JOIN(I)
300      CONTINUE
3000  CONTINUE
C--
C  That's it.
C--
      STOP
      END
