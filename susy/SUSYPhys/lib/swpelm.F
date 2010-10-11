C----------------------------------------------------------------------C
      SUBROUTINE SWPELM( ARRAY1, IORDER, ARRAY2, IDIM )
C----------------------------------------------------------------------C
C*  Swap ARRAY1(*) in decreasing order.
C*  Input  ARRAY1(*), IDIM
C*  Output IORDER(*), ARRAY2(*)
C----------------------------------------------------------------------C
      REAL*4    ARRAY1(IDIM), ARRAY2(IDIM)
      INTEGER*4 IORDER(IDIM)
C--------------------------------------
C  Reset arrays first.
      DO 1 I = 1, IDIM
        IORDER(I) = 0
1       ARRAY2(I) = 0
 
      DO 10 I1 = 1, IDIM
       IORDER(I1) = 1
       DO 20 I2 = 1, IDIM
        IF( ARRAY1(I1) .LT. ARRAY1(I2) ) IORDER(I1) = IORDER(I1)+1
20     CONTINUE
       ARRAY2(IORDER(I1)) = ARRAY1(I1)
10    CONTINUE
 
      DO 30 I = 2, IDIM
        IYES = 0
        DO 40 J = 1, IDIM
         IF( I .EQ. IORDER(J) ) IYES = 1
40      CONTINUE
      IF( IYES .EQ. 0 ) ARRAY2(I) = ARRAY2(I-1)
30    CONTINUE
 
C    Check IORDER(*)
      DO 50 I = 1, IDIM
       IF( IORDER(I) .LE. 0 .OR. IORDER(I) .GT. IDIM ) THEN
         WRITE(6,*) ' Error at SWPELM: Order array (IORDER)',
     &              ' may be wrong.'
         STOP
       ENDIF
50    CONTINUE
 
      RETURN
      END
