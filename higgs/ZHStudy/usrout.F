      SUBROUTINE USROUT

      CHARACTER*128 BSHSTF
C--
C  Set logical units and open the histogram file.
C--
      LUHIST = 24
      BSHSTF = 'bases.tdr'
      OPEN(LUHIST, FILE=BSHSTF)
C--
C  Output histograms.
C-- 
      DO 100 IH = 1, 12
         CALL XHSAVE( LUHIST,  IH )
100   CONTINUE       
C--
C  That's it.
C-- 
      RETURN
      END
