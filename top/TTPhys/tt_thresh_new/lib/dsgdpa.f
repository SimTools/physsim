CC**********************************************************************
C*
C*===========================================================
C* Subroutine DSGDPA(MODE,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
C*                                              EE,QR, DSG)
C*=====================================================---===
C*
C* (Purpose)
C*    Calculates Green_q function to the ttbar threshold shape
C*    determined from the EW amplitude for e+e- --> bW+ bbarW-.
C* (Inputs)
C*    MODE   : (I*4) : mode.
C*    ALFS   : (R*8) : alpha_s.
C*    ALF    : (R*8) : alpha.
C*    S2W    : (R*8) : sin^2(theta_W).
C*    AMZ    : (R*8) : m_Z.
C*    AMW    : (R*8) : m_W.
C*    AMB    : (R*8) : m_b.
C*    AMT    : (R*8) : m_t.
C*    VTB2   : (R*8) : !V_tb!^2
C*    EE     : (R*8) : sqrt(s)-2m_t.
C*    QR     : (R*8) : !p_top!.
C* (Output)
C*    DSG    : (R*8) : ds/dp.
C* (Update Record)
C*   93/04/19  K.Fujii        Original version.
C*
CC**********************************************************************
 
      SUBROUTINE DSGDPA(MODE,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                                              EE,QR, DSG)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE
      REAL*8       ALFS, ALF, S2W, AMZ, AMW, AMB, AMT, VTB2,
     .             EE, QR, DSG
C--
      PARAMETER    ( NE = 100, NQ = 500 )
      INTEGER*4    NBN(3)
      REAL   *8    ERGN(0:3), DERGN(3), DQRGN(0:NE,3)
      REAL   *8    DSGL, DSGH, DSDAT(0:NQ,0:NE,3)
C>>>  For m_t = 150 GeV.
      DATA NBN     /   5, 35, 25 /
      DATA ERGN    / -12.D0, -7.D0, 0.D0, 5.D0 /
C>>>
      SAVE
C--
C  Statement function.
C--
C      BETA(X1,X2) = SQRT( MAX( 1 - 2*(X1+X2) + (X1-X2)**2, 0.D0 ) )
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      MODEP = MODE
      IF ( MODEP.LE.2 ) THEN
         PRINT *, ' ****** DSGDPA initialization '
         PRINT *, '   DSGDPA now tabulates ds/dp for'
         PRINT *, '     ALFS   = ', ALFS
         PRINT *, '     AMT    = ', AMT
         PRINT *, '     VTB2   = ', VTB2
         DO 100 IRGN = 1, 3
            ME  = NBN(IRGN)
            EMN = ERGN(IRGN-1)
            EMX = ERGN(IRGN)
            DE  = (EMX-EMN)/ME
            DERGN(IRGN) = DE
            DO 10 IE = 0, ME
               E     = EMN + DE*IE
               QMAX  = SQRT( (3*AMT+E-AMW)*(3*AMT+E+AMW)
     .                 *(AMT+E-AMW)*(AMT+E+AMW) )/2/(2*AMT+E)
               DQ    = QMAX/NQ
               DQRGN(IE,IRGN) = DQ
               DO 1 IQ = 0, NQ
                  Q  = IQ*DQ
                  CALL DSGDP0(MODEP,ALFS,ALF,S2W,AMZ,AMW,AMB,
     .                        AMT,VTB2,  E,Q, DSG )
                  DSDAT(IQ,IE,IRGN) = DSG
                  MODEP = 4
1              CONTINUE
               MODEP = 3
10          CONTINUE
100      CONTINUE
         PRINT *, 'ds/dp tabulation ended.'
      ENDIF
C--
C  Calculate ds/dp.
C--
      E    = EE
      IF ( E.LE.ERGN(0) ) THEN
         DSG = 0
         RETURN
      ELSE IF ( E.GT.ERGN(3) ) THEN
         PRINT *, '>>>>> Error in DSGDPA >>>>>>>>'
         PRINT *, '    E = ', EE, ' is too large!'
         STOP
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
            DSGL = DSDAT(IQL,IE ,IRGN)
     .              + (DSDAT(IQL+1,IE ,IRGN)-DSDAT(IQL,IE ,IRGN))*FL
         ELSE
            DSGL = DSDAT(NQ ,IE ,IRGN)
         ENDIF
C--
         DQH  = DQRGN(IE+1,IRGN)
         IQH  = QR/DQH
         IF ( IQH.LT.NQ ) THEN
            FH   = (QR-IQH*DQH)/DQH
            DSGH = DSDAT(IQH,IE+1,IRGN)
     .             + (DSDAT(IQH+1,IE+1,IRGN)-DSDAT(IQH,IE+1,IRGN))*FH
         ELSE
            DSGH = DSDAT(NQ ,IE+1,IRGN)
         ENDIF
         DSG = DSGL + (DSGH-DSGL)*F
      ELSE
         DQ  = DQRGN(ME,IRGN)
         IQ  = QR/DQ
         F   = (QR-IQ*DQ)/DQ
         DSG = DSDAT(IQ,ME,IRGN)
     .          + (DSDAT(IQ+1,ME,IRGN)-DSDAT(IQ,ME,IRGN))*F
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
