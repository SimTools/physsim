      SUBROUTINE GQTABL(ALFS,AMT,AMH,BTH,RS,QR,GQ)

      IMPLICIT    REAL*8 ( A-H, O-Z )
      REAL*8      ALFS, AMT, AMH, BTH, RS, QR
      COMPLEX*16  GQ
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
      DATA NBN     /   23, 35, 300 /
C     DATA NBN     /   46, 35, 600 /
C     DATA NBN     /  230, 140, 600 /
C<<< 20090817
      REAL   *8    ERGN(0:3), DERGN(3), DQRGN(0:NE,3)
      DATA ERGN    / -30.D0, -7.D0, 0.D0, 60.D0 /
      COMPLEX*16   GQL, GQH, GQDAT(0:NQ,0:NE,3)
      DATA NCALL   /0/
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
C--
C  Calculate Green function.
C--
      E    = RSH - 2*AMST
      IF ( E.LE.ERGN(0) ) THEN
         GQ = 0
         RETURN
      ELSE IF ( E.GT.ERGN(3) ) THEN
C         PRINT *, '>>>>> Error in QCDMOD >>>>>>>>'
C         PRINT *, '   RS = ', RS, ' is too large!'
C         STOP

C     no QCD correction for large E
         GQ = 0
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
#if 0
         print *, ' ie = ', ie, ' irgn = ', irgn, ' iql = ', iql,
     .            ' iqh = ', iqh, ' dql = ', dql, ' dqh = ', dqh
         print *, ' gq(iql,ie) = ', gqdat(iql,ie,irgn)
         print *, ' gq(iqh,ie) = ', gqdat(iqh,ie,irgn)
#endif
      ELSE
         DQ  = DQRGN(ME,IRGN)
         IQ  = QR/DQ
         F   = (QR-IQ*DQ)/DQ
         GQ  = GQDAT(IQ,ME,IRGN)
     .          + (GQDAT(IQ+1,ME,IRGN)-GQDAT(IQ,ME,IRGN))*F
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
