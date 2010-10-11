 
      SUBROUTINE SGSMSM(RS,POLE,AM,CF,Q,T3,SG,R)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      PARAMETER  ( IxPRC = 2 )
      REAL   *8  RS, POLE, AM, CF, Q, T3, SG, R
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xPI    = ACOS(-1.D0)
         x2PI   = 2*xPI
         x4PI   = 4*xPI
         xSIN2W = 0.23
         xCOS2W = 1 - xSIN2W
         xSINW  = SQRT(xSIN2W)
         xCOSW  = SQRT(xCOS2W)
         xCOTW  = xCOSW/xSINW
         xALF   = 1.D0/128D0
         xGV2PB = 3.893796623D8
         AMZ    = 91.17D0
         AMZ2   = AMZ*AMZ
C--
         C0     = SQRT(x4PI*xALF)
         CV     = C0/(xSINW*xCOSW)*(-0.25D0+xSIN2W)
         CA     = C0/(xSINW*xCOSW)*( 0.25D0 )
         FAC    =  xGV2PB/(48*xPI)
      ENDIF
C--
C  Set coupling constants, etc.
C--
         CG     = C0*Q
         CZ     = C0/(xSINW*xCOSW)*(T3-Q*xSIN2W)
C--
C  Set independent variables.
C--
      EBM  = RS/2
      BT   = (EBM-AM)*(EBM+AM)
      IF ( BT.LE.0.D0 ) THEN
         SG = 0
         RETURN
      ELSE
         BT   = SQRT(BT)/EBM
      ENDIF
      S    = RS*RS
      A    = -CG*C0/S + CV*CZ/(S-AMZ2)
      B    =            CA*CZ/(S-AMZ2)
C--
C  Set independent variables.
C--
      PL  = (1-POLE)/2
      PR  = (1+POLE)/2
      SG  = FAC*((A-B)**2*PL+(A+B)**2*PR)*S*BT**3*CF
      SG0 = xGV2PB*x4PI/3*xALF**2/S
      R   = SG/SG0
C--
C  That's it.
C--
      RETURN
      END
