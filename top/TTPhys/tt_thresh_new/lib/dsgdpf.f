 
      SUBROUTINE DSGDPF(MODE,ALFS,ALF,S2W,AMZ,AMW,AMB,
     .                  AMT, VTB2,  E,Q,DSGEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      INTEGER*4    MODE
      REAL*8       ALFS, ALF, S2W, AMZ, AMW, AMB, AMT, VTB2,
     .             E, Q, DSGEFF
      REAL*4       XRS
      DATA         NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xPI    = ACOS(-1.D0)
         AME    = 0.511D-3
         ALF0   = 1./137.
         FACT   = 1 + 2*ALF0/xPI * ( xPI**2/6 - 1.D0/4 )
      ENDIF
C--
C  Set bin size for numerical integration.
C--
      IF ( AMT.LT.130.D0 ) THEN
         NBM    = 200
      ELSE
         NBM    = 100
      ENDIF
      DBM    = 1.D0/NBM
C--
C  Gaussian beam width, beamstrahlung, and bremsstrahlng.
C--
      MODEP  = MODE
      RS     = 2*AMT + E
      DSGEFF = 0
      DO 10 IBM = 0, NBM
         IF ( IBM.EQ.0 .OR. IBM.EQ.NBM ) THEN
            WGT = 1.D0/3
         ELSE IF ( MOD(IBM,2).EQ.0 ) THEN
            WGT = 2.D0/3
         ELSE
            WGT = 4.D0/3
         ENDIF
         CALL RSHDIS(REAL(DBM*IBM),1,XRS)
         EH = RS*XRS - 2*AMT
C--
         CALL DSGDPA(MODEP,ALFS,ALF,S2W,AMZ,AMW,AMB,AMT,VTB2,
     .                                            EH,Q, DSG)
         DSGEFF = DSGEFF + DSG*WGT*DBM
         MODEP  = 3
10    CONTINUE
C--
C  Vertex correction.
C--
      DSGEFF = FACT*DSGEFF
C--
C  That's it.
C--
      RETURN
      END
