      SUBROUTINE IDLPTN(LNxDT,INDAT,ID)
 
      INTEGER*4   LNxDT, ID, INDAT(LNxDT)
      INTEGER*4   IDATA(20)
      REAL   *4   RDATA(20)
      EQUIVALENCE (IDATA,RDATA)
C
C==================<<  Entry Point  >>================================
C
C--
C  Copy input to local array.
C--
      CALL UVCOPY(20,INDAT(1),IDATA(1))
C--
C  Look at track nature.
C--
      ICHG = IDATA( 9)
      IDET = IDATA(10)
      NCDC = IDATA(11)
      NEMC = IDATA(12)
      NHDC = IDATA(13)
      IF ( ICHG.EQ.0 ) THEN
         ID = 0
         RETURN
      ELSE
         IF ( NEMC.GT.0 .AND. NHDC.EQ.0 ) ID = -ICHG*11
         IF ( NEMC.EQ.0 .AND. NHDC.EQ.0 ) ID = -ICHG*13
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
