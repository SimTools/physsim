C***********************************************************************
C*
C*  -------------------------------------------------====
C*  Subroutine UDSRCH( NREC, LREC, KEY, ARAY, XDATA, NPNT )
C*  -------------------------------------------------====
C*(Function)
C*   Make binary search.    Double presision version.
C*   Find NPNT which satisfies
C*        ARAY(KEY,NPNT) < or = XDATA < ARAY(KEY,NPNT+1)
C*   ARAY should be sorted from smaller data to larger data.
C*
C*(Input)
C*   NREC  ; # of records
C*   LREC  ; Record size  (ARAY(LREC,NREC))
C*   KEY   ; KEy element pointer in a record.
C*   ARAY  ; Target data aray
C*   XDATA ; Reference data.
C*(Output)
C*   NPNT  ; Pointer which satisfies the conditions.
C*           NPNT = 0 if XDATA < ARAY(KEY,1)
C*           NPNT = NREC if XDATA =ARAY(KEY,NREC) or
C*                                > ARAY(KEY,NREC)
C*(Author)
C*   A. Miyamoto   16-Oct-1990  Original version.
C*
C***********************************************************************
C
      SUBROUTINE UDSRCH(NREC,LREC,KEY,ARAY,XDATA,NPNT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8   ARAY(LREC,NREC)
C
C ====< Entry Point > ==================================================
C
C ------------------------------------
C (1) Omit Extream case first.
C ------------------------------------
C
      IF( XDATA .LT. ARAY(KEY,1) ) THEN
        NPNT = 0
        RETURN
      ELSEIF( XDATA .GE. ARAY(KEY,NREC) ) THEN
        NPNT = NREC
        RETURN
      ENDIF
C
C
C ------------------------------------
C (2) Binary search
C ------------------------------------
C
C
      I1 = 1
      I2 = NREC
200   CONTINUE
      IMDL = (I2+I1)/2
      IF( XDATA.GE.ARAY(KEY,IMDL) ) THEN
          I1 = IMDL
      ELSE
          I2 = IMDL
      ENDIF
C
      IF( I2-I1 .GT. 1 ) GO TO 200
      NPNT = I1
      RETURN
      END
