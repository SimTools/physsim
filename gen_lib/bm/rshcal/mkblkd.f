C*  91/08/31  K.Fujii                  Original version.
C*                                     This version makes a block
C*                                     data for RSHDIS from data
C*                                     created by RSHCAL.
C*
      IMPLICIT    REAL*8   ( A-H, O-Z )
      PARAMETER ( MXxPT = 500 )
      REAL*4      RDATA(2,0:MXxPT)
      DATA LIU    /10/
      DATA LOU    /20/
C--
C  Read in sqrt(s) distribution.
C--
      NP = -1
1     READ(LIU,*,END=10) X, Y
      NP = NP + 1
      RDATA(1,NP) = X
      RDATA(2,NP) = Y
      GO TO 1
C--
C  End of data.
C--
10    PRINT *, ' NP  = ', NP
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
         IF ( IS.GT.NS ) LL = NL/2
         IF ( MOD(LL,10).NE.0 ) LL = LL + 1
         KST = MOD(LST+1,2)
         DO 20 L = LF, LL
            IP = I + 2*(L-1)
            IF ( L.NE.LL ) THEN
               WRITE(LOU,'(5X,''.'',4(D15.8,'',''))')
     .                               RDATA(2,IP), RDATA(1,IP),
     .                               RDATA(2,IP+1), RDATA(1,IP+1)
            ELSE
               IF ( KST.EQ.0 ) THEN
                  WRITE(LOU,'(5X,''.'',3(D15.8,'',''),D15.8,''/'')')
     .                                  RDATA(2,IP), RDATA(1,IP),
     .                                  RDATA(2,IP+1), RDATA(1,IP+1)
               ELSE
                  WRITE(LOU,'(5X,''.'',D15.8,'','',D15.8,''/'')')
     .                                  RDATA(2,IP), RDATA(1,IP)
               ENDIF
            ENDIF
20       CONTINUE
200   CONTINUE
C--
C  That's it.
C--
      STOP
      END
