C        IMPLICIT REAL*8 ( A-H, O-Z )
C        ALF = 1/128.D0
C        S2W = 0.230D0
C        AMW = 80
C        AMB = 5
C        VTB = 1
C        PRINT *, ' Type in m_t.'
C  1     READ(5,*,END=999) AMT
C        IF ( AMT.LT.0.D0 ) GO TO 999
C        CALL GTGAMT(ALF,S2W,AMW,AMB,AMT,VTB,GAMT)
C        PRINT *, ' m_t = ', AMT, ' GeV  Gamma_t = ', GAMT, ' GeV'
C        GO TO 1
C  999   STOP
C        END
 
      SUBROUTINE  GTGAMT(ALF,S2W,AMW,AMB,AMT,VTB,GAMT)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ALF, S2W, AMW, AMB, AMT, VTB, GAMT
C
C========< Entry Point >================================================
C
C--
C  Calculate on-shell width.
C--
      XMBR   = (AMB/AMT)**2
      XMWR   = (AMW/AMT)**2
      PFACT  = DSQRT(1.0D0 + (XMWR-XMBR)**2 - 2.0D0*(XMWR + XMBR))
      BFACT  = (1.0D0- XMBR)**2 + (1.0D0 + XMBR)*XMWR - 2.0D0*XMWR**2
      GAMT   = ((ALF/S2W)/16)*(AMT**3/AMW**2)*(VTB**2)*PFACT*BFACT
C--
C  That's it.
C--
      RETURN
      END
