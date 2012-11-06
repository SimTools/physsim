CC**********************************************************************
C*
C*===================================================----==
C* Subroutine QCDMOD(ALFS,AMT,GMT,AMH,BTH,RS,Q12,Q22,FQCD)
C*===================================================----==
C*
C* (Purpose)
C*    Calculates QCD modification to the ttbar threshold shape
C*    determined from the EW amplitude for e+e- --> bW+ bbarW-.
C* (Inputs)
C*    ALFS   : (R*4) : alpha_s.
C*    AMT    : (R*4) : top mass.
C*    GMT    : (R*4) : top width.
C*    AMH    : (R*4) : Higgs mass.
C*    BTH    : (R*4) : top's Yukawa coupling normalized by SM.
C*    RS     : (R*4) : sqrt(s).
C*    Q12    : (R*4) : virtual top mass squared
C*    Q22    : (R*4) : virtual anti-top mass squared.
C* (Output)
C*    FQCD   : (C*8) : (Q1,Q2)-dependent QCD modification to
C*                   : the ttbar threshold calculated without QCD.
C* (Update Record)
C*   92/10/27  K.Fujii        Original version.
C*                            This version allows variable sqrt(s)
C*                            in the range 2*mt-12 to 2*mt+5,
C*                            though parameters are not allowed to
C*                            change in the job.
C*   92/10/27  K.Fujii        QCD hard gluon correction added.
C*   07/11/08  A.Ishikawa     Allowed sqrt(s) is changed to 
C*                            2*mt-20 to 2*mt+40GeV for ttH.
C*                            Set FQCD=1 for sqrt(s)>2*mt+40GeV.
C*   09/08/13  K.Ikematsu     Extended QCD correction to 2*mt+60GeV.
C*   09/08/20  K.Fujii        Modified for better matching with the
C*                            analytic results 
C*
CC**********************************************************************
 
      SUBROUTINE QCDMOD(ALFS,AMT,GMT,AMH,BTH,RS,Q12,Q22,FQCD)
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
      REAL*4      ALFS, AMT, GMT, AMH, BTH, RS, Q12, Q22
      COMPLEX*8   FQCD
C--
#if 1
C     PARAMETER    ( NE = 300, NQ = 500 )
      PARAMETER    ( NE = 600, NQ = 500 )
#else
C     PARAMETER    ( NE = 300, NQ = 5000 )
      PARAMETER    ( NE = 600, NQ = 5000 )
#endif
      INTEGER*4    NBN(3)
C>>> 20090817
C     DATA NBN     /   13, 35, 300 /
      DATA NBN     /   46, 35, 600 /
C<<< 20090817
      REAL   *8    ERGN(0:3), DERGN(3), DQRGN(0:NE,3)
      DATA ERGN    / -30.D0, -7.D0, 0.D0, 60.D0 /
      COMPLEX*16   GQ, GQL, GQH, GQDAT(0:NQ,0:NE,3)
      DATA NCALL   /0/
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT( MAX( 1 - 2*(X1+X2) + (X1-X2)**2, 0.D0 ) )
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         PRINT *, ' ****** QCDMOD initialization '
         PRINT *, '    Can''t change parameters! '
         NCALL  = 1
         PI     = ACOS(-1.D0)
         ALPS   = ALFS
         ALP    = 1/128.D0
         SN2W   = 0.23D0
         AMSZ   = 91.17D0
         GAMZ   = 2.5D0
         AMSW   = 80.D0
         AMSB   =  5.D0
         AMST   =  AMT
         GAMT   =  GMT
         AMSH   =  AMH
         BETH   =  BTH
         VFF    =  1
         QCDVTX =  ABS( 1 - 8*ALFS/3/PI )
C--
         PRINT *, 'QCDMOD now tabulates GQ for'
         PRINT *, '  ALFS   = ', ALPS
         PRINT *, '  AMT    = ', AMST
         PRINT *, '  VTB2   = ', VFF**2
         PRINT *, '  AMH    = ', AMH
         PRINT *, '  BTH    = ', BTH
         PRINT *, '  GAMT   = ', GAMT
         MODE  =  0
         DO 100 IRGN = 1, 3
            ME  = NBN(IRGN)
            EMN = ERGN(IRGN-1)
            EMX = ERGN(IRGN)
            DE  = (EMX-EMN)/ME
            DERGN(IRGN) = DE
            DO 10 IE = 0, ME
               E     = EMN + DE*IE
               QMAX  = SQRT( (3*AMST+E-AMSW)*(3*AMST+E+AMSW)
     .                 *(AMST+E-AMSW)*(AMST+E+AMSW) )/2/(2*AMST+E)
               DQ    = QMAX/NQ
               DQRGN(IE,IRGN) = DQ
               RSH   = E + 2*AMST
               DO 1 IQ = 0, NQ
                  Q  = IQ*DQ
                  CALL GREENQ(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                       AMST, VFF, AMSH,BETH,   RSH,Q,GQ)
                  GQDAT(IQ,IE,IRGN) = GQ
                  MODE = 3
1              CONTINUE
               MODE = 2
10          CONTINUE
100      CONTINUE
         PRINT *, 'GQ tabulation ended.'
      ENDIF
C--
C  Calculate the magnitude of the 3-momentum of t and tbar quarks.
C--
      RSH  = RS
      S    = RSH*RSH
      QR   = (RSH/2)*BETA(Q12/S,Q22/S)
C--
C  Calculate Green function.
C--
      E    = RSH - 2*AMST
      IF ( E.LE.ERGN(0) ) THEN
         FQCD = 1
         RETURN
      ELSE IF ( E.GT.ERGN(3) ) THEN
C         PRINT *, '>>>>> Error in QCDMOD >>>>>>>>'
C         PRINT *, '   RS = ', RS, ' is too large!'
C         STOP

C     no QCD correction for large E
         FQCD = 1
         RETURN
      ELSE IF ( E.LE.ERGN(1) ) THEN
         IRGN = 1
      ELSE IF ( E.LE.ERGN(2) ) THEN
         IRGN = 2
      ELSE
         IRGN = 3
      ENDIF
      ME  = NBN(IRGN)
      EMN = ERGN(IRGN-1)
      DE  = DERGN(IRGN)
C--
      IE  = (E-EMN)/DE
      IF ( IE.LT.ME ) THEN
         EMN = ERGN(IRGN-1)
         DE  = DERGN(IRGN)
         F   = (E-EMN-IE*DE)/DE
C--
         DQL  = DQRGN(IE  ,IRGN)
         IQL  = QR/DQL
         IF ( IQL.LT.NQ ) THEN
            FL   = (QR-IQL*DQL)/DQL
            GQL  = GQDAT(IQL,IE  ,IRGN)
     .              + (GQDAT(IQL+1,IE  ,IRGN)-GQDAT(IQL,IE  ,IRGN))*FL
         ELSE
            GQL  = GQDAT(NQ ,IE  ,IRGN)
         ENDIF
C--
         DQH  = DQRGN(IE+1,IRGN)
         IQH  = QR/DQH
         IF ( IQH.LT.NQ ) THEN
            FH   = (QR-IQH*DQH)/DQH
            GQH  = GQDAT(IQH,IE+1,IRGN)
     .              + (GQDAT(IQH+1,IE+1,IRGN)-GQDAT(IQH,IE+1,IRGN))*FH
         ELSE
            GQH  = GQDAT(NQ ,IE+1,IRGN)
         ENDIF
         GQ  = GQL + (GQH-GQL)*F
      ELSE
         DQ  = DQRGN(ME,IRGN)
         IQ  = QR/DQ
         IF( IQ.LT.NQ ) THEN 
           F= (QR-IQ*DQ)/DQ
           GQ  = GQDAT(IQ,ME,IRGN)
     .          + (GQDAT(IQ+1,ME,IRGN)-GQDAT(IQ,ME,IRGN))*F
         ELSE
           GQ = GQDAT(NQ,ME,IRGN)
         ENDIF 
      ENDIF
C--
C  Calculate Green function.
C--
#if 0
      FQCD = ( RSH - 2*AMST - QR*QR/AMST + (0.D0,1.D0)*GAMT )*GQ
#else
      RET  = AMST/SQRT(AMST**2+QR**2)
      FQCD = ( (Q12 + Q22)/(2*AMST) - AMST  + (0.D0,1.D0)*GAMT*RET )*GQ
#endif
      W    = (AMSW/AMST)**2
      RW   = AMST**3*(1+2*W)*(1-W)**2
      Q1   = SQRT(Q12)
      X1   = (AMSW/Q1)**2
      IF (X1.LT.0.9) THEN
         R1   = RW/(Q1**3*(1+2*X1)*(1-X1)**2)
      ELSE
         R1   = 1
      ENDIF
      Q2   = SQRT(Q22)
      X2   = (AMSW/Q2)**2
      IF (X2.LT.0.9) THEN
         R2   = RW/(Q2**3*(1+2*X2)*(1-X2)**2)
      ELSE
         R2   = 1
      ENDIF    
C--
      FQCD = FQCD*SQRT(R1*R2)
      FQCD = QCDVTX*FQCD
C--
C  That's it.
C--
      RETURN
      END
