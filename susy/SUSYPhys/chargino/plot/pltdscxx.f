CINCLUDE DSGCXX
CINCLUDE XCMASS
C*
C*  Differential cross section for e+e- --> X+ X-.
C*    RS = sqrt(S)
C*
C*  5/16/91  K.Fujii       This version can handle X+1 X-2 case.
C*                         Tree level differential cross section.
C*
      IMPLICIT    REAL*8  ( A-H, O-Z )
      REAL   *8   GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2),
     .            AMSN, GMSN, AMZ, GMZ, AMX(2), GMX(2), AMXC(2)
      REAL   *8   SG(0:3)
      REAL   *8   UL(2,2), UR(2,2)
C--
      CHARACTER*20 PATRN(1:4)
      DATA PATRN  / '0.2 0.0',
     .              '0.2 0.05 0.01 0.05',
     .              '0.2 0.05',
     .              '0.01 0.05 0.01' /
      CHARACTER*9 JOIN(1:4)
      DATA JOIN   / 'SOLID  ',
     .              'PATTERNED',
     .              'PATTERNED',
     .              'PATTERNED'/
C--
C     SUSY(1,*,*) = mu
C         (2,*,*) = m_2
C         (3,*,*) = tan(beta)
C         (4,*,*) = m_snu
C--
      REAL   *8   SUSY(4,3,8)
      DATA SUSY /
     .  400.,180.,+50.,248.,  250.,250.,+2., 547.,  250.,500.,+2., 547.,     
C    .   21., 10.,+1.,13.37,  250.,250.,+2., 547.,  250.,500.,+2., 547.,
C    .  500.,250.,+2., 547.,  250.,250.,+2., 547.,  250.,500.,+2., 547.,     
C    .  500.,200.,-2., 500.,  200.,200.,-2., 500.,  200.,500.,-2., 500.,
     .  500.,200.,+2., 500.,  200.,200.,+2., 500.,  200.,500.,+2., 500.,
     .  500.,200.,-8., 500.,  200.,200.,-8., 500.,  200.,500.,-8., 500.,
     .  500.,200.,+8., 500.,  200.,200.,+8., 500.,  200.,500.,+8., 500.,
     .  500.,200.,-2.,1000.,  200.,200.,-2.,1000.,  200.,500.,-2.,1000.,
     .  500.,200.,+2.,1000.,  200.,200.,+2.,1000.,  200.,500.,+2.,1000.,
     .  500.,200.,-8.,1000.,  200.,200.,-8.,1000.,  200.,500.,-8.,1000.,
     .  500.,200.,+8.,1000.,  200.,200.,+8.,1000.,  200.,500.,+8.,1000./
      DATA LOU / 20 /
C
C========< Entry Point >================================================
C
C--
C  Differential cross section.
C--
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET WINDOW X 2 10 Y 2 8'')')
C     WRITE(LOU,'(''SET LIMITS X -1. 1. Y 0. 0.025'')')
      WRITE(LOU,'(''SET LIMITS X -1. 1. Y 0. 0.100'')')
      WRITE(LOU,'(''(SET LIMITS X -1. 1. Y 0.001 10.'')')
      WRITE(LOU,'(''(SET SCALE Y LOGARITHMIC'')')
      WRITE(LOU,'(''SET LABELS SIZE 2.5'')')
      WRITE(LOU,'(''SET TITLE  SIZE 4.0'')')
      WRITE(LOU,'(''TITLE   6.5 7.5 SIZE 2.5'',
     .            '' '''' e2+3e2-3 --> C013S242+3C013S242-3'')')
      WRITE(LOU,'(''CASE                   '',
     .            '' ''''  X X X X     GUVVMVVX XGUVVMVVX X'')')
      WRITE(LOU,'(''TITLE 5.7 1.2 ''''cos Q'')')
      WRITE(LOU,'(''CASE          ''''    G'')')
      WRITE(LOU,'(''TITLE 0.2 3.7 ANGLE 90.'',
     .                            '' ''''dS/dW (pb)'')')
      WRITE(LOU,'(''CASE          '''' G  F     '')')
      WRITE(LOU,'(''SET ORDER X Y'')')
C--
C  Define constants.
C--
      xPI       = ACOS(-1.)
      x2PI      = 2*xPI
      x4PI      = 2*x2PI
C--
C  Standard Model parameters.
C--
      xSIN2W    = 0.23
      xCOS2W    = 1 - xSIN2W
      xSINW     = SQRT(xSIN2W)
      xCOSW     = SQRT(xCOS2W)
      xALF      = 1.D0/128D0
      AMZ       = 91.1
      AMW       = 80.0
      GMZ       = 2.5
C--
C  Angular region for plot.
C     CSMN = cos(theta)_min
C     CSMX = cos(theta)_max
C--
      CSMN = -1
      CSMX = -CSMN
      NCS  = 200
      DCS  = (CSMX-CSMN)/NCS
      PHI  = 0.
C--
C  Select reaction type and parameter set.
C     IXM  = chi- : (1,2) = (chi1,chi2)
C     IXP  = chi+ : (1,2) = (chi1,chi2)
C     ISET = susy parameter set 
C--
      IXM       = 1
      IXP       = 1
      ISET      = 1
C--
C  Set sqrt(s) and polarization.
C--
      RS        = 400.
      POL       = -0.9
C--
C  Couplings of electron to Z and gamma.
C--
      C0        = SQRT(x4PI*xALF)
      CW        = C0/xSINW
      CZ        = CW/xCOSW
      QE        = -1
      T3LE      = -0.5D0
      T3RE      =  0
C--
      GAL(1)    = C0*QE
      GAL(2)    = C0*QE
      GZL(1)    = CZ*(T3LE-QE*xSIN2W)
      GZL(2)    = CZ*(T3RE-QE*xSIN2W)
C--
      WRITE(LOU,'(''TITLE 7.2 7.0 SIZE 2.5'',
     .            '' ''''2s0O = '',F4.1,'' TeV'')')   RS/1000
      WRITE(LOU,'(''CASE          ''''M UD      '')')
      WRITE(LOU,'(''TITLE 2.5 7.5 SIZE 2.'',
     .            '' ''''m0N013S241   = '',F5.0,'' GeV'')')
     .            SUSY(4,1,ISET)
      WRITE(LOU,'(''CASE '''' XGUVVMVVX   '')')
      WRITE(LOU,'(''TITLE 2.5 7.2 SIZE 2.'',
     .                '' ''''tanB  = '',F5.0)') SUSY(3,1,ISET)
      WRITE(LOU,'(''CASE ''''   G '')')
C--
      WRITE(LOU,'(''TITLE 2.5 6.9 SIZE 2.'',
     .            '' ''''(M,M021) = ('',F5.0,'','',F5.0,''):'')')
     .            SUSY(1,1,ISET), SUSY(2,1,ISET)
      WRITE(LOU,'(''CASE '''' G  X X   '')')
      WRITE(LOU,'(''SET ORDER X Y'')')
C      WRITE(LOU,'(''5.7 6.9;6.5,6.9;JOIN TEXT SOLID'')')
      WRITE(LOU,'(''TITLE 5.7 6.9 SIZE 2. ''''solid'')')
      WRITE(LOU,'(''TITLE 2.5 6.6 SIZE 2.'',
     .            '' ''''       = ('',F5.0,'','',F5.0,''):'')')
     .            SUSY(1,2,ISET), SUSY(2,2,ISET)
      WRITE(LOU,'(''SET ORDER X Y'')')
C      WRITE(LOU,'(''5.7 6.6;6.5,6.6;JOIN TEXT DOTDASH'')')
      WRITE(LOU,'(''TITLE 5.7 6.6 SIZE 2. ''''dotdash'')')
      WRITE(LOU,'(''TITLE 2.5 6.3 SIZE 2.'',
     .            '' ''''       = ('',F5.0,'','',F5.0,''):'')')
     .            SUSY(1,3,ISET), SUSY(2,3,ISET)
      WRITE(LOU,'(''SET ORDER X Y'')')
C      WRITE(LOU,'(''5.7 6.3;6.5,6.3;JOIN TEXT DASH'')')
      WRITE(LOU,'(''TITLE 5.7 6.3 SIZE 2. ''''dash'')')
C--
C  Loop over SUSY parameters.
C--
      DO 100 ICASE = 1, 3
C--
C  Set SUSY parameters.
C--
         AMU       = SUSY(1,ICASE,ISET)
         AM2       = SUSY(2,ICASE,ISET)
         TNB       = SUSY(3,ICASE,ISET)
         BT        = MOD(ATAN(TNB)+xPI,xPI)
         AMSN      = SUSY(4,ICASE,ISET)
         GMX(1)    = 0.
         GMX(2)    = 0.
         GMSN      = 0.
C--
C  Calculate masses and mixings.
C--
         CALL XCMASS(AM2,AMU,BT,AMW,AMXC,FIL,FIR,EPSR)
         AMX(1) = AMXC(IXM)
         AMX(2) = AMXC(IXP)
C--
C  Coupling of chargino to Z and gamma.
C     U_i = ! cos(phi_i)  -sin(phi_i) !
C           ! sin(phi_i)   cos(phi_i) !
C  where i = (L,R). U_i transforms weak eigen states(wino,higgsino)
C  to mass eigen states(charginos):
C   ! X_- !   = !  cos(phi_i)  sin(phi_i) ! * ! W^- !
C   ! X_+ !_i   ! -sin(phi_i)  cos(phi_i) !   ! H^- !_i
C--
         CSL       = COS(FIL)
         SNL       = SIN(FIL)
         CSR       = COS(FIR)
         SNR       = SIN(FIR)
C--
         UL(1,1)   =  CSL
         UL(1,2)   =  SNL
         UL(2,1)   = -SNL
         UL(2,2)   =  CSL
C--
         UR(1,1)   =  CSR
         UR(1,2)   =  SNR
         UR(2,1)   = -SNR*EPSR
         UR(2,2)   =  CSR*EPSR
C--
         QX        = -1
         T3LX      =   UL(IXM,1)*UL(IXP,1)*(-1)
     .               + UL(IXM,2)*UL(IXP,2)*(-1./2)
         T3RX      =   UR(IXM,1)*UR(IXP,1)*(-1)
     .               + UR(IXM,2)*UR(IXP,2)*(-1./2)
C--
         GAX(1)    = C0*QX*(UL(IXM,1)*UL(IXP,1)+UL(IXM,2)*UL(IXP,2))
         GAX(2)    = C0*QX*(UR(IXM,1)*UR(IXP,1)+UR(IXM,2)*UR(IXP,2))
         GZX(1)    = CZ*(T3LX-QX*xSIN2W)
         GZX(2)    = CZ*(T3RX-QX*xSIN2W)
         GSNX(1,1) = CW*UR(IXM,1)
         GSNX(2,1) = 0
         GSNX(1,2) = 0
         GSNX(2,2) = CW*UR(IXP,1)
C--
         PRINT *, ' GAX       = ', GAX
         PRINT *, ' GZX       = ', GZX
         PRINT *, ' GSNX(*,1) = ', GSNX(1,1), GSNX(2,1)
         PRINT *, ' GSNX(*,2) = ', GSNX(1,2), GSNX(2,2)
C--
         SGSM = 0
         WRITE(LOU,'(''(                                    '')')
         WRITE(LOU,'(''(  **********************************'')')
         WRITE(LOU,'(''(    AMU   = '',F10.3)') SUSY(1,ICASE,ISET)
         WRITE(LOU,'(''(    AM2   = '',F10.3)') SUSY(2,ICASE,ISET)
         WRITE(LOU,'(''(    TANB  = '',F10.3)') SUSY(3,ICASE,ISET)
         WRITE(LOU,'(''(    AMSN  = '',F10.3)') SUSY(4,ICASE,ISET)
         WRITE(LOU,'(''(   ---------------------------------'')')
         WRITE(LOU,'(''(    PHIL  = '',F10.3)') FIL
         WRITE(LOU,'(''(    PHIR  = '',F10.3)') FIR
         WRITE(LOU,'(''(    EPS   = '',F10.3)') EPSR
         WRITE(LOU,'(''(    AMM   = '',F10.3)') AMX(1)
         WRITE(LOU,'(''(    AMP   = '',F10.3)') AMX(2)
         WRITE(LOU,'(''(    AMX1  = '',F10.3)') AMXC(1)
         WRITE(LOU,'(''(    AMX2  = '',F10.3)') AMXC(2)
         WRITE(LOU,'(''(  **********************************'')')
         WRITE(LOU,'(''(                                    '')')
C--
         DO 10   ICS = 0, NCS
            IF ( ICS.EQ.0 .OR. ICS.EQ.NCS ) THEN
               WGT = 1./3
            ELSE IF ( MOD(ICS,2).EQ.0 ) THEN
               WGT = 2./3
            ELSE
               WGT = 4./3
            ENDIF
            CS = CSMN + DCS*ICS
            CALL DSGCXX(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .                  AMX, GMX ,RS, POL, CS, PHI, SG)
            WRITE(LOU,'(2E15.5)') CS, SG(0)
            SGSM = SGSM + SG(0)*WGT
10       CONTINUE
         WRITE(LOU,'(''SET PATTERN '',A20)') PATRN(ICASE)
         WRITE(LOU,'(''JOIN '',A9)') JOIN(ICASE)
         WRITE(LOU,'(''( sqrt<s> ='',F8.1,
     .                  '' sigma = '',E15.5,'' pb'')')
     .                  RS, SGSM*x2PI*DCS
100   CONTINUE
C--
      STOP
      END
