C*  8/30/91  K.Fujii       Convolute bremsstrahlung and beam energy
C*                         spread and output reduced sqrt(s) dist.
C*                         This version ignores beamstrahlung.
C*
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   SDATA(2,0:1000)
      REAL*8   RDATA(4,0:1000)
      DATA LOU /20/ LOU2 /10/
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C     ITYP = 1 : flat
C          = 2 : Gaussian
C          = 3 : Realistic
C--
C>>>
      ITYP  = 3
      ENSG  = 1/2.D0
      FG    = 1.2D0
      SGEB  = 1.D-6
C     SGEB  = 0.0005
C     SGEB  = 0.0010
C     SGEB  = 0.0020
C     SGEB  = 0.0035
C     SGEB  = 0.0050
C     SGEB  = 0.0070
C>>>
      xALF  = 1/137.D0
      xPI   = ACOS(-1.D0)
      AME   = 0.511D-3
      AMT   = 150
      ROOTS = 2*AMT
C--
C  First create sqrt(s) distribution with beam energy spread.
C--
      NFR  = 500
      FRMN = 1 - SGEB
      FRMX = 1 + SGEB
      DFR  = ( FRMX - FRMN )/NFR
      NZ    = NFR
      DO 1000 IFR = 0, NFR
         FR   = FRMN + IFR*DFR
         PRB  = 0
         Z1MN = MAX( ( 1 + SGEB - FR**2/(1-SGEB) )/2/SGEB, 0.D0 )
         Z1MX = MIN( ( 1 + SGEB - FR**2/(1+SGEB) )/2/SGEB, 1.D0 )
         DZ   = (Z1MX-Z1MN)/NZ
         DO 100 IZ1 = 0, NZ
            Z1 = Z1MN + IZ1*DZ
            FM = 1 + SGEB - 2*SGEB*Z1
            FP = FR**2/FM
            Z2 = ( 1 + SGEB - FP )/2/SGEB
            IF ( ITYP.EQ.1 ) THEN
               PR1  = 1
               PR2  = 1
            ELSE IF ( ITYP.EQ.2 ) THEN
               PR1  = EXP(-((2*Z1-1)/ENSG)**2/2 )
               PR2  = EXP(-((2*Z2-1)/ENSG)**2/2 )
            ELSE
               PR1  = 1 - EXP(-((2*Z1-1)/ENSG/2)**2/2 )/FG
               PR2  = 1 - EXP(-((2*Z2-1)/ENSG/2)**2/2 )/FG
            ENDIF
            IF ( IZ1.EQ.0 .OR. IZ1.EQ.NZ ) THEN
               W1 = 1
            ELSE IF ( MOD(IZ1,2).EQ.0 ) THEN
               W1 = 2
            ELSE
               W1 = 4
            ENDIF
            PRB = PRB + PR1*PR2*DZ*W1/3
100      CONTINUE
         SDATA(1,IFR) = FR
         SDATA(2,IFR) = PRB
C        PRINT *, ' FR = ', FR, ' PRB = ', PRB
1000  CONTINUE
C--
C  Loop over fraction.
C--
      EPS   = 1.D-8
      FRCMX = 1 + SGEB
      FRCMN = 1.D-2
      NFRAC  = 100
      FRACMN = - LOG( 1 + SGEB - FRCMN + EPS )
      FRACMX = - LOG( 1 + SGEB - FRCMX + EPS )
      DFRAC  = ( FRACMX - FRACMN )/NFRAC
C--
      SUM = 0
      DO 2000 IFRAC = 1, NFRAC
         FRAC  = FRACMN + (IFRAC  )*DFRAC
         FRACL = FRACMN + (IFRAC-1)*DFRAC
         FRC   = 1 + SGEB + EPS - EXP(-FRAC)
         FRCL  = 1 + SGEB + EPS - EXP(-FRACL)
         DFRC  = FRC - FRCL
C--
C  Start integration.
C     Z1 : Energy spread of e- beam.
C      2 : Energy spread of e+ beam.
C--
         PRB   = 0
         NFRP  = NFR
         FRPMX = FRMX
         FRPMN = MAX(FRCL,FRMN)
         DFRP  = (FRPMX-FRPMN)/NFRP
C>>>
         DO 200 IFRP = 0, NFRP
C        DO 200 IFRP = 0, 0
C           FRPMN = 1
C>>>
            FRP  = FRPMN + IFRP*DFRP
            IFR  = (FRP-FRMN)/DFR
            IFR  = MIN(IFR,NFR-1)
            PR   = SDATA(2,IFR) +
     .             (SDATA(2,IFR+1)-SDATA(2,IFR))*(FRP-SDATA(1,IFR))/DFR
            IF ( IFRP.EQ.0 .OR. IFRP.EQ.NFRP ) THEN
               W1 = 1
            ELSE IF ( MOD(IFRP,2).EQ.0 ) THEN
               W1 = 2
            ELSE
               W1 = 4
            ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
            RS   = ROOTS*FRP
            BTE  = (2*xALF/xPI)*(2*LOG(MAX(RS,20*AME)/AME)-1)
            ZG   = ( 1 - FRC/FRP )*( 1 + FRC/FRP )
            ZGL  = MAX(0.D0,ZG)
            ZGH  = ( 1 - FRCL/FRP )*( 1 + FRCL/FRP )
            ADD  = ( (ZGH**BTE-ZGL**BTE)*(1+3*BTE/4)
     .             - BTE*(ZGH-ZGL)*(1-(ZGH+ZGL)/4) )/DFRC
            ADD  = ADD*PR*DFRP*(W1/3)
            PRB  = PRB + ADD
200      CONTINUE
C--
C  End of integration.
C--
         SUM = SUM + PRB*DFRC
         RDATA(1,IFRAC) = FRAC
         RDATA(2,IFRAC) = SUM
         RDATA(3,IFRAC) = FRC
         RDATA(4,IFRAC) = PRB
2000  CONTINUE
C--
C  Prepare topdraw file.
C--
C     WRITE(LOU,*) 'NEW FRAME'
C     WRITE(LOU,*) 'SET SCALE Y LOG'
C     WRITE(LOU,*) 'SET ORDER X Y'
      DO 300 IFRAC = 0, NFRAC
         RDATA(2,IFRAC) = RDATA(2,IFRAC)/SUM
         WRITE(LOU2,*) RDATA(1,IFRAC), RDATA(2,IFRAC)
300   CONTINUE
C     WRITE(LOU,*) 'JOIN'
C--
C  Prepare topdraw file.
C--
      WRITE(LOU,*) 'NEW FRAME'
      WRITE(LOU,*) 'SET SCALE Y LOG'
      WRITE(LOU,*) 'SET ORDER X Y'
      DO 400 IFRAC = 0, NFRAC
         RDATA(4,IFRAC) = RDATA(4,IFRAC)/SUM
         WRITE(LOU,*) RDATA(3,IFRAC), RDATA(4,IFRAC)
C        WRITE(LOU,*) RDATA(1,IFRAC), RDATA(4,IFRAC)
400   CONTINUE
      WRITE(LOU,*) 'JOIN'
C--
C  That's it.
C--
      STOP
      END
