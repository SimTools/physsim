C
C ********************************************************
C
C      SUBROUTINE MHIGGS(AMA, TANBET, FM3, AT, AB, XMU,
C     >            SFM3M, ALF, S2W, AMZ, 
C     >            ALPH,AMHL,AMHH,AMHC )
C
C(FUNCTION)
C   CALCULATE THE MASS OF HIGGS BOSON IN ON-LOOP LEVEL
C   Call EWINIT prior to call this routine.  EWINIT is attached below
C
C(INPUT)
C   AMA  : PSEUDO SCALAR HIGGS MASS.
C   TANBET  : TAN(BETA) = V2/V1
C   FM3     :(MTAU,MTOP,MBOTTOM)
C   AT      : TRILINEAR SOFT SUPERSYMMETRY-BREAKING PARAMETER FOR TOP
C   AB      : TRILINEAR SOFT SUPERSYMMETRY-BREAKING PARAMETER FOR BOTTOM
C   XMU     : HIGGSINO MASS
C (OUTPUT)
C   AMHL   : LIGHT HIGGS MASS
C   AMHH   : HEAVY HIGGS MASS 
C   AMHC   : CHARGED HIGGS MASS ALA TREE LEVEL MASS RELATION 
C   ALPH    : HIGGS MIXING ANGLE.
C(AUTHOR)
C   A. MIYAMOTO  8-APR-1991  ORIGINAL VERSION.
C(REFERENCE)
C   J.ELLIS, ET AL., PL 262B(1991)477
C
C **********************************************************************
C23456
      SUBROUTINE MHIGGS(AMA, TANBET, FM3, AT, AB, XMU,
     >            SFM3M, ALF, S2W, AMZ, 
     >            ALPH,AMHL,AMHH,AMHC )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  PI, S2W, AMZ
      REAL*8 FM3(3),SFM3M(7)
      GFUNC(X1,X2) = 2.D0 - (X1+X2)/(X1-X2)*LOG(X1/X2)
C
C =====< ENTRY POINT >==================================================
C
      AMT=FM3(2)
      AMB=FM3(3)
      BETA = ATAN(TANBET)
      COSB = COS(BETA)
      SINB = SIN(BETA)
      COS2B = COSB*COSB - SINB*SINB
      COSBSQ = COSB*COSB
      SINBSQ = SINB*SINB
      PI = ACOS(-1.D0)
      GWEAK = SQRT(ALF*4*PI/S2W)
      AMZ2 = AMZ*AMZ
      AMW2=AMZ2*(1-S2W)

C
C (1) MASS OF SCALAR PARTICLES
C
C .. FOR SCALOR TOP
      STMS1=SFM3M(4)**2
      STMS2=SFM3M(5)**2
C
C .. FOR SCALOR BOTTOM
      SBMS1=SFM3M(6)**2
      SBMS2=SFM3M(7)**2     
 
C
C  ONE LOOP CORRECTION
C
      AMB4 = AMB**4
      AMT4 = AMT**4
      BBTERM  = AB*(AB+XMU*TANBET)/(SBMS1 - SBMS2)
      BTTERM  = XMU*(AT+XMU/TANBET)/(STMS1-STMS2)
      TTTERM  = AT*(AT+XMU/TANBET)/(STMS1 - STMS2)
      TBTERM  = XMU*(AB+XMU*TANBET)/(SBMS1-SBMS2)
      GBS     = GFUNC(SBMS1, SBMS2)
      GTS     = GFUNC(STMS1, STMS2)
      DLOGB   = LOG(SBMS1/SBMS2)
      DLOGT   = LOG(STMS1/STMS2)
     
       D11  = AMB4/COSBSQ*( LOG(SBMS1*SBMS2/AMB4) + 2*BBTERM*LOGB )
C       write(*,*) d11
       d11  =AMB4/COSBSQ*BBTERM**2*GBS
C       write(*,*) d11
       d11 = AMT4/SINBSQ*BTTERM**2*GTS
C       write(*,*) d11
      D11  = AMB4/COSBSQ*( LOG(SBMS1*SBMS2/AMB4) + 2*BBTERM*LOGB )
     >     + AMB4/COSBSQ*BBTERM**2*GBS
     >     + AMT4/SINBSQ*BTTERM**2*GTS
C       write(*,*) amt4, sinbsq, btterm, gts, d11
      D22  = AMT4/SINBSQ*( LOG(STMS1*STMS2/AMT4) + 2*TTTERM*DLOGT )
     >       + AMT4/SINBSQ*TTTERM**2*GTS
     >       + AMB4/COSBSQ*TBTERM**2*GBS
 
      D12  = AMT4/SINBSQ*BTTERM*( DLOGT + TTTERM*GTS )
     >       + AMB4/COSBSQ*TBTERM*( DLOGB + BBTERM*GBS )
C
C  MASS MATRIX.
C
      DELTA = AMA**2*2*SINB*COSB
      C1    = AMZ2*2*SINB*COSB
      C3    = 3*GWEAK**2/(8*PI**2*AMW2)
C      C3 = 0
C      write(*,*) c3
      A11   = C1/TANBET + TANBET*DELTA + C3*D11
      A12   = -1*C1 - DELTA + C3*D12
      A22   = C1*TANBET + DELTA/TANBET + C3*D22
      DSQ   = (A11+A22)**2 - 4*(A11*A22-A12**2)
 
      IF( DSQ.LT.0.D0 ) THEN
            PRINT *,'IN HIGGS..D <0.  SOMETHING MUS BE WRONG.'
        RETURN
      ENDIF
C      write(*,*) a11,a22, a11+a22,sqrt(dsq)
      AMH1S    = 0.25D0*( (A11+A22) - SQRT(DSQ) )
      AMH2S    = 0.25D0*( (A11+A22) + SQRT(DSQ) )
C      write(*,*) amh1s,amh2s
      AMHL = SQRT(AMH1S)
      AMHH = SQRT(AMH2S)
      TANALP = (AMH1S+AMH2S)*A12/( AMH2S*A11 - AMH1S*A22 )
      ALPH = ATAN(TANALP)
      AMHC = SQRT(AMA**2+AMW2)
 
      RETURN
      END
