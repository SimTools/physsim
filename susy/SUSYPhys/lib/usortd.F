C*********************************************************************
C*
C*    ----------------------------------------======
C*    SUBROUTINE USORTD(NREC, LREC, KEY, ARY , IPNT)
C*    ----------------------------------------======
C*
C*(Function)
C*    Quick sorter of an real*4 array.
C*
C*(Input)
C*    NREC    :  Number of records.
C*    LREC    :  Record size.
C*    KEY     :  Key element pointer in a record.
C*    ARY     :  Real    array to be sortted.
C*
C*(Output)
C*    IPNT    :  Pointer containing the results.
C*
C*(Author)
C*
C*    Hirofumi Fujii    Jan 12, 1983
C*    A. Miyamoto         6-Sep-1985, modified the arguments.
C*
C*********************************************************************
C*
      SUBROUTINE USORTD(NREC, LREC, KEY, ARY, IPNT)
C
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4   IPNT(NREC)
      REAL*8      ARY(LREC, NREC)
      PARAMETER  (ISTKSZ = 32 ,
     >            ILEFT  = 1,
     >            IRIGHT = 2 )
      INTEGER     ISTACK(2, ISTKSZ)
C
C*
C*(1)  Initialize.
C*
      DO 100 I = 1, NREC
         IPNT(I) = I
100   CONTINUE
      ISP   = 1
      ISTACK(ILEFT ,ISP) = 1
      ISTACK(IRIGHT,ISP) = NREC
C*
C*(2)  Get the most left and most right position from the stack.
C*
200   CONTINUE
        NLEFT  = ISTACK(ILEFT, ISP)
        NRIGHT = ISTACK(IRIGHT,ISP)
        ISP    = ISP - 1
C*
C*(3)  Determine the range to be sortted.
C*
300     CONTINUE
          IL   = NLEFT
          IR   = NRIGHT
          IM   = (IL+IR)/2
C>>>
C           VALKEY = ARY(KEY, IPNT(IM))
          VALKEY = ABS(ARY(KEY, IPNT(IM)))
C>>>
C*
C*(4)   Swap the data with respect to the VALKEY.
C*
400       CONTINUE
C>>>
C 410         IF(ARY(KEY,IPNT(IL)).GE.VALKEY) GO TO 420
410         IF(ABS(ARY(KEY,IPNT(IL))).GE.VALKEY) GO TO 420
C>>>
              IL = IL + 1
              IF(IL.GT.NREC)
     >        PRINT *,'Warning in USORTD  IL greater than NREC...'
 
              GO TO 410
C>>>
C  420        IF(ARY(KEY,IPNT(IR)).LE.VALKEY) GO TO 430
420         IF(ABS(ARY(KEY,IPNT(IR))).LE.VALKEY) GO TO 430
C>>>
              IR = IR - 1
              IF(IR.LE.0)
     >        PRINT *,'Warning in USORTD  IR less than 0....'
              GO TO 420
C
430         CONTINUE
C
C    Swap data.
C
            IF(IL.LE.IR) THEN
              ISAVE = IPNT(IL)
              IPNT(IL) = IPNT(IR)
              IPNT(IR) = ISAVE
              IL = IL + 1
              IR = IR - 1
            END IF
 
            IF(IL.LE.IR) GO TO 400
 
C*
C*(5)  Find out the next region to be sorted.
C*
            IF(IR-NLEFT.LT.NRIGHT-IL) THEN
              IF(IL.LT.NRIGHT) THEN
                ISP = ISP + 1
                ISTACK(ILEFT ,ISP) = IL
                ISTACK(IRIGHT,ISP) = NRIGHT
              ENDIF
              NRIGHT = IR
            ELSE
              IF(NLEFT.LT.IR) THEN
                ISP = ISP + 1
                ISTACK(ILEFT ,ISP) = NLEFT
                ISTACK(IRIGHT,ISP) = IR
              ENDIF
              NLEFT = IL
            ENDIF
 
          IF(NLEFT.LT.NRIGHT) GO TO 300
 
        IF(ISP.GT.0) GO TO 200
 
      RETURN
      END
