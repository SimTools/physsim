C*       IMPLICIT  REAL*8 ( A-H, O-Z )
C*       PARAMETER    ( MXxCNT = 3, MXxLM = 60, MXxAM = 6 )
C*       REAL*8       RDATA(3,0:MXxLM,0:MXxAM), CHICNT(MXxCNT)
C*       DATA LOU     / 20 /
C*
C* C--
C*       CHI2MN = 6.02
C*       CHI1S  = CHI2MN + 1.00
C*       CHI68  = CHI2MN + 2.28
C*       CHI90  = CHI2MN + 4.61
C*       CHICNT(1) = CHI1S
C*       CHICNT(2) = CHI68
C*       CHICNT(3) = CHI90
C* C--
C*       DO 10 IL = 0, MXxLM
C*          DO 1 IM = 0, MXxAM
C*             READ(10,*) (RDATA(I,IL,IM),I=1,3)
C* 1        CONTINUE
C* 10    CONTINUE
C* C--
C*       CALL UCONT2(LOU,MXxCNT,CHICNT,MXxLM,MXxAM,RDATA)
C* C--
C*       STOP
C*       END
C*
C*
C*
C*
      SUBROUTINE UCONT2(LOU,MXxCNT,CHICNT,MXxLM,MXxAM,RDATA)
 
      IMPLICIT     REAL*4 ( A-H, O-Z )
      INTEGER*4    LOU, MXxCNT, MXxLM, MXxAM
      REAL*8       RDATA(3,0:MXxLM,0:MXxAM), CHICNT(MXxCNT)
C--
      PARAMETER    ( MXxPNT = 100, MXxPAR = 4 )
      REAL*4       XDAT(0:MXxPNT), YDAT(0:MXxPNT), DYDAT(0:MXxPNT),
     .             A(MXxPAR), DA(MXxPAR), XYDAT(3,MXxPNT,2), CHIF
C
C========< Entry Point >================================================
C
C--
C  Set frame.
C--
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET LIMITS X '',2E12.5,'' Y '',2E12.5)')
     .           RDATA(1,0,0), RDATA(1,MXxLM,MXxAM),
     .           RDATA(2,0,0), RDATA(2,MXxLM,MXxAM)
      WRITE(LOU,'(''SET ORDER X Y DUMMY'')')
C--
C  First scan through X variable values and decide shift and scale
C  parameters.
C--
      XMN =  1.E20
      XMX = -1.E20
      DO 100 ILM = 0, MXxLM
         DO 10 IAM = 0, MXxAM
            X = RDATA(1,ILM,IAM)
            IF ( X.LT.XMN ) XMN = X
            IF ( X.GT.XMX ) XMX = X
10       CONTINUE
100   CONTINUE
      XSHIFT = (XMX+XMN)/2
      XSCALE = 1/(XMX-XMN)
C--
C  Then modify X variable.
C--
      DO 200 ILM = 0, MXxLM
         DO 20 IAM = 0, MXxAM
            RDATA(1,ILM,IAM) = (RDATA(1,ILM,IAM)-XSHIFT)*XSCALE
20       CONTINUE
200   CONTINUE
C--
C  Loop over contours.
C--
      DO 9000 ICNT = 1, MXxCNT
         CHI2  = CHICNT(ICNT)
         NPNTL = 0
         NPNTH = 0
         IAMS  = MXxAM/2
         DO 700 ILM = 0, MXxLM
            CHI2S = 1.E20
            ALAM  = RDATA(2,ILM,0)
            DO 30 IAM = 0, MXxAM
               XDAT(IAM)  = RDATA(1,ILM,IAM)
               YDAT(IAM)  = RDATA(3,ILM,IAM)
               DYDAT(IAM) = 1
C>>>
C     PRINT *,   XDAT(IAM), ALAM, YDAT(IAM)
C>>>
               IF ( CHI2S.GT.RDATA(3,ILM,IAM) ) THEN
                  IAMS  = IAM
                  CHI2S = RDATA(3,ILM,IAM)
               ENDIF
30          CONTINUE
C--
C  Now check if any solution there.
C--
            IF ( IAMS.EQ.0 ) THEN
               IF ( RDATA(3,ILM,IAMS).GT.CHI2 )  GO TO 700
               IAMMN = IAMS
               IAMMX = MIN( IAMS+4, MXxAM )
            ELSE IF ( IAMS.EQ.MXxAM ) THEN
               IF ( RDATA(3,ILM,IAMS).GT.CHI2 )  GO TO 700
               IAMMN = MAX( IAMS-4, 0 )
               IAMMX = IAMS
            ELSE
               IF ( IAMS.LT.2 ) THEN
                  IAMMN = 0
                  IAMMX = MIN( IAMMN+4, MXxAM )
               ELSE
                  IAMMX = MIN( IAMS+2, MXxAM )
                  IAMMN = MAX( IAMMX-4, 0 )
               ENDIF
            ENDIF
            NAM = IAMMX - IAMMN + 1
C--
C  Parabola fit the chi**2 curve for this slice.
C--
            CALL POLFIT(XDAT(IAMMN),YDAT(IAMMN),DYDAT(IAMMN),
     .                  NAM,3,1,A,DA,CHIF)
            A0 = A(1)
            A1 = A(2)
            A2 = A(3)
C--
C  Redo POLFIT.
C--
            CALL POLFIT(XDAT(IAMMN),YDAT(IAMMN),DYDAT(IAMMN),
     .                  NAM,MXxPAR,1,A,DA,CHIF)
C--
C  No solution. Go to next slice.
C--
            IF ( RDATA(3,ILM,IAMS).GT.CHI2 ) THEN
               XMN = -A1/2/A2
               YMN = A0 - A2*XMN**2
               IF ( YMN.GT.CHI2 )                GO TO 700
C--
C  Find approximate lower solution.
C--
               B0  = A0 - CHI2
               B1  = A1
               B2  = A2
               SQRTD =  B1**2 - 4*B2*B0
               AMT0L = ( -B1 - SQRTD )/2/B2
               AMT0H = ( -B1 + SQRTD )/2/B2
            ELSE
               DO 40 IAM = 0, MXxAM-1
                  CHI2L = RDATA(3,ILM,IAM  )
                  CHI2H = RDATA(3,ILM,IAM+1)
                  IF ( CHI2.LE.CHI2L .AND. CHI2.GE.CHI2H ) THEN
                     AMTL  = RDATA(1,ILM,IAM  )
                     AMTH  = RDATA(1,ILM,IAM+1)
                     AMT0L = AMTL
     .                      + (CHI2-CHI2L)*(AMTH-AMTL)/(CHI2H-CHI2L)
                                                 GO TO 45
                  ENDIF
40             CONTINUE
C--
45             DO 50 IAM = MXxAM, 1, -1
                  CHI2L = RDATA(3,ILM,IAM-1)
                  CHI2H = RDATA(3,ILM,IAM  )
                  IF ( CHI2.GE.CHI2L .AND. CHI2.LE.CHI2H ) THEN
                     AMTL  = RDATA(1,ILM,IAM-1)
                     AMTH  = RDATA(1,ILM,IAM  )
                     AMT0H  = AMTL
     .                      + (CHI2-CHI2L)*(AMTH-AMTL)/(CHI2H-CHI2L)
                                                 GO TO 55
                  ENDIF
50             CONTINUE
            ENDIF
C--
C  Now look for solution by Newton's method.
C--
55          AMT0 = AMT0L
            CALL USOLVE(MXxPAR,A,AMT0,CHI2,AMT,CHI2S,IRT)
            NPNTL = NPNTL + 1
            IF ( IRT.NE.0 ) THEN
               XYDAT(1,NPNTL,1) = AMT0
               XYDAT(2,NPNTL,1) = ALAM
               XYDAT(3,NPNTL,1) = -2
            ELSE
               XYDAT(1,NPNTL,1) = AMT
               XYDAT(2,NPNTL,1) = ALAM
               XYDAT(3,NPNTL,1) = CHI2S
            ENDIF
C--
            AMT0 = AMT0H
            CALL USOLVE(MXxPAR,A,AMT0,CHI2,AMT,CHI2S,IRT)
            NPNTH = NPNTH + 1
            IF ( IRT.NE.0 ) THEN
               XYDAT(1,NPNTH,2) = AMT0
               XYDAT(2,NPNTH,2) = ALAM
               XYDAT(3,NPNTH,2) = -2
            ELSE
               XYDAT(1,NPNTH,2) = AMT
               XYDAT(2,NPNTH,2) = ALAM
               XYDAT(3,NPNTH,2) = CHI2S
            ENDIF
700      CONTINUE
C--
C  Now draw contours.
C--
         DO 80 IP = 1, NPNTL
C           IF ( IP.GT.2 .AND. IP.LT.NPNTL-1 .AND. MOD(IP,6).NE.0 )
C    .       GOTO 80
            XYDAT(1,IP,1) = XYDAT(1,IP,1)/XSCALE + XSHIFT
            WRITE(LOU,'(3E20.8)') (XYDAT(I,IP,1),I=1,3)
80       CONTINUE
         DO 90 IP = NPNTH, 1, -1
C           IF ( IP.GT.2 .AND. IP.LT.NPNTH-1 .AND. MOD(IP,6).NE.0 )
C    .       GOTO 90
            XYDAT(1,IP,2) = XYDAT(1,IP,2)/XSCALE + XSHIFT
            WRITE(LOU,'(3E20.8)') (XYDAT(I,IP,2),I=1,3)
90       CONTINUE
         WRITE(LOU,'(3E20.8)') (XYDAT(I,1,1),I=1,3)
         WRITE(LOU,'(''SET SYMBOL 9O SIZE .2; PLOT'')')
         WRITE(LOU,'(''JOIN'')')
9000  CONTINUE
C--
C  That's it.
C--
      RETURN
      END
