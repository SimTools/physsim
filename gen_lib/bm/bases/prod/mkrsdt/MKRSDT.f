C* For BASES
      IMPLICIT    REAL*8   ( A-H, O-Z )
      PARAMETER ( MXxPT = 200 )
      REAL*4      RDATA(0:3,0:MXxPT)
      DATA LIU    /5/
      DATA LOU    /6/
C--
C  Read in sqrt(s) distribution.
C--
      NP = 0
      RDATA(1,0)  = 0
1     READ(LIU,*,END=10) X, Y, DY
      NP = NP + 1
      RDATA(0,NP) = X - RDATA(1,NP-1)
      RDATA(1,NP) = X
      RDATA(2,NP) = Y
      RDATA(3,NP) = DY
      GO TO 1
C--
C  End of data.
C--
10    PRINT *, ' NP  = ', NP
      RDATA(0,1) = RDATA(1,2) - RDATA(1,1)
      RDATA(0,0) = RDATA(0,1)
C--
C  Decide where bin width changed.
C--
      DO 11 IP = 1, NP
         DXL = RDATA(0,IP-1)
         DXR = RDATA(0,IP)
         IF ( ABS(DXL-DXR).GT.1.E-5 )            GO TO 12
11    CONTINUE
12    NPB = IP - 1
      PRINT *, ' NPB = ', NPB
C--
C  Integrate Y from 0 to X.
C--
      SUM = 0
      DX = RDATA(1,NPB) - RDATA(1,NPB-1)
      DO 15 IP = 1, NPB
         SUM = SUM + RDATA(2,IP)*DX
         RDATA(0,IP) = SUM
15    CONTINUE
C--
      DX = RDATA(1,NPB+2) - RDATA(1,NPB+1)
      DO 16 IP = NPB+1, NP
         SUM = SUM + RDATA(2,IP)*DX
         RDATA(0,IP) = SUM
16    CONTINUE
      PRINT *, ' SUM = ', SUM
C--
      DO 17 IP = 1, NP
         RDATA(0,IP) = RDATA(0,IP)/SUM
17    CONTINUE
C--
C  Prepare block data.
C--
      NS    = NP/20
      NL    = NP - 20*NS
      ISLST = 1
      IF ( NL.EQ.0 ) ISLST = 0
      DO 200 IS = 1, NS+ISLST
         IF ( IS.LE.NS ) THEN
            LST = 19
         ELSE
            LST = NL - 1
         ENDIF
         I = 1 + 20*(IS-1)
         J = I + LST
         WRITE(LOU,'(''      DATA ((BMDATA(I,J),I=1,2),J='',I3,
     .               '','',I3,'') /'')') I, J
         LF = 1
         LL = 10
         IF ( IS.GT.NS ) LL = NL/2 + MOD(NL,2)
         KST = MOD(LST+1,2)
         DO 20 L = LF, LL
            IP = I + 2*(L-1)
            IF ( L.NE.LL ) THEN
               WRITE(LOU,'(5X,''.'',4(D15.8,'',''))')
     .                               RDATA(0,IP), RDATA(1,IP),
     .                               RDATA(0,IP+1), RDATA(1,IP+1)
            ELSE
               IF ( KST.EQ.0 ) THEN
                  WRITE(LOU,'(5X,''.'',3(D15.8,'',''),D15.8,''/'')')
     .                                  RDATA(0,IP), RDATA(1,IP),
     .                                  RDATA(0,IP+1), RDATA(1,IP+1)
               ELSE
                  WRITE(LOU,'(5X,''.'',D15.8,'','',D15.8,''/'')')
     .                                  RDATA(0,IP), RDATA(1,IP)
               ENDIF
            ENDIF
20       CONTINUE
200   CONTINUE
C--
C  That's it.
C--
      STOP
      END
