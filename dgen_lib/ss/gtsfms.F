C* (Update Record)
C*    92/06/09  K.Fujii      Original version.
C*                           This version does not take into account
C*                           L-R mixings.
C*
 
      SUBROUTINE GTSFMS(AM0,AM2,AMU,TNB,ALF,S2W,AMZ,SFM)
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
      REAL*8      AM0, AM2, AMU, TNB, ALF, S2W, SFM(7)
      REAL*8      ALFT(3), B(3), F(3), C(2,7), AMSF2(7)
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants at Q**2 = MZ**2
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         AMU   = AMU
         x2PI  = 2*ACOS(-1.D0)
         x4PI  = 2*x2PI
         NGN   = 3
C--MSSM
         B(1)  = (10/3.D0)*NGN + 1
         B(2)  = 2*NGN - 6 + 1
         B(3)  = 2*NGN - 9
      ENDIF
C--
C  Constants at Q**2 = MZ**2
C--
      AMZCS = AMZ*AMZ*(TNB-1)*(TNB+1)/(TNB*TNB+1)
C--
C  MX.
C--
      TEE   = (x4PI/ALF)*(3-8*S2W)/(3*B(1)-5*B(2))
      AMX   = AMZ*EXP(TEE/2)
      ALFWI = S2W/ALF
      ALFXI = ALFWI - (B(2)/x4PI)*TEE
      ALFSI = ALFXI + (B(3)/x4PI)*TEE
C--
C     PRINT *, ' 1/ALF(MZ) = ', 1/ALF
C     PRINT *, ' SIN2W     = ', S2W
C     PRINT *, ' AMZ       = ', AMZ, ' GeV'
C     PRINT *, ' NGN       = ', NGN
C     PRINT *, ' AMX       = ', AMX
C     PRINT *, ' ALFS      = ', 1/ALFSI
C     PRINT *, ' --- '
C--
C  Alpha_tilde.
C--
      ALFT(3) = 1/ALFXI/x4PI
      ALFT(2) = ALFT(3)
      ALFT(1) = 3/5.D0*ALFT(2)
C--
C  f_i.
C--
      DO 10 I = 1, 3
         F(I) = 1/B(I)/ALFT(I) * ( 1 - 1/(1+B(I)*ALFT(I)*TEE)**2 )
10    CONTINUE
C--
C  M_i**2 = M_0**2 + C(1,i)*M**2 + C(2,i)*cos(2*beta)*M_Z**2.
C--
      C(1,1) = 2*( 4/3.D0*ALFT(3)*F(3) + 3/4.D0*ALFT(2)*F(2)
     .                               + 1/36.D0*ALFT(1)*F(1) )
      C(1,2) = C(1,1)
      C(1,3) = 2*( 4/3.D0*ALFT(3)*F(3) + 4/9.D0*ALFT(1)*F(1) )
      C(1,4) = 2*( 4/3.D0*ALFT(3)*F(3) + 1/9.D0*ALFT(1)*F(1) )
      C(1,5) = 2*( 3/4.D0*ALFT(2)*F(2) + 1/4.D0*ALFT(1)*F(1) )
      C(1,6) = C(1,5)
      C(1,7) = 2*ALFT(1)*F(1)
C--
      C1QSM  = C(1,1) + C(1,2) + C(1,3) + C(1,4)
      C1LSM  = C(1,5) + C(1,6) + C(1,7)
C--
      C(2,1) = -0.5D0 + 2/3.D0*S2W
      C(2,2) =  0.5D0 - 1/3.D0*S2W
      C(2,3) =      - 2/3.D0*S2W
      C(2,4) =      + 1/3.D0*S2W
      C(2,5) =  0.5D0 -      S2W
      C(2,6) = -0.5D0
      C(2,7) =      +      S2W
C--
C  Loop over sferminons. AMSF2 is mass**2.
C     AMSF2(I) = (SUL,SDL,SUR,SDR,SEL,SNL,SER)
C--
      AM   = ALFWI/ALFXI*AM2
      AM02 = AM0**2
      DO 20 I = 1, 7
         AMSF2(I) = AM02 + AM*AM*C(1,I) + C(2,I)*AMZCS
20    CONTINUE
C--
      SFM(1) = SQRT(AMSF2(6))
      SFM(2) = SQRT(AMSF2(5))
      SFM(3) = SQRT(AMSF2(7))
      SFM(4) = SQRT(AMSF2(1))
      SFM(5) = SQRT(AMSF2(3))
      SFM(6) = SQRT(AMSF2(2))
      SFM(7) = SQRT(AMSF2(4))
C--
C  That's it.
C--
      RETURN
      END
