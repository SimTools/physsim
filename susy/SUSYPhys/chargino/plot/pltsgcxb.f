CX*INCLUDE RSHDIS
CINCLUDE RSHDIS0
CINCLUDE UDSRCH
C--
CINCLUDE SGCXXA
C--
C*INCLUDE XCMASS
CINCLUDE COUPLSUB
CINCLUDE INOMIX
CINCLUDE SWPELM
CINCLUDE INSUSY
CINCLUDE DCY2BODY
CINCLUDE GTSFMS
CINCLUDE SF2BR
CINCLUDE HGSMIX
C--
C*
C*
C*  Total cross section for e+e- --> X+ X-.
C*    RS = sqrt(S)
C*
C*  5/16/91  K.Fujii       This version can handle X+1 X-2 case.
C*			   This version calculates cross section
C*			   with ISR and beam effects.
C*  1/11/95  K.Fujii	   Brought from FACOM to jlcux.       
C*			   The beam effects can be switched off
C*                         by setting BEAMLIB = gbm_0_lib.a in
C*			   the make file.
C*
      IMPLICIT    REAL*4  ( A-H, O-Z )
      REAL   *4  AM0, AMU, AM2, TNB, ALF, ALFS, S2W, ZM, FM(3),
     .           SFM(7), SWM(2), SZM(4), GMSF(7), GMSW(2), GMSZ(4)
      COMPLEX*8  GNCW(2,4,2), GCCZ(2,2,2), GNNZ(2,4,4)
      COMPLEX*8  GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     .           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*8  GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     .           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
C--
      REAL   *4  GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2),
     .           AMSN, GMSN, AMZ, GMZ, AMX(2), GMX(2)
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
C         (4,*,*) = m_0
C--
      REAL   *4   SUSY(4,3,8)
      DATA SUSY /
C    .  400.,150.,+2., 400.,  150.,150.,+2., 400.,  150.,400.,+2., 400.,
CSTD     .  400.,250.,+2.,  70.,  250.,250.,+2.,  70.,  250.,400.,+2.,  70.,
C    .  400.,250.,+3.,  70.,  250.,250.,+3.,  70.,  250.,400.,+2.,  70.,
     .  400.,250.,+3., 500.,  250.,250.,+3.,  70.,  250.,400.,+2.,  70.,
C    .  400.,250.,+3.,  70.,  250.,250.,+3.,  70.,  250.,400.,+3.,  70.,
CKatoy     .  400.,180.,+50.,200.,  250.,250.,+3.,  70.,  250.,400.,+3.,  70.,
C    .  400.,250.,+10., 70.,  250.,250.,+10., 70.,  250.,400.,+10., 70.,
C    .  400.,250.,+2.,  20.,  400.,250.,+2.,  70.,  400.,250.,+2., 120.,
C Point 3 .  263.,83.7,-2., 200.,  400.,250.,+2.,  70.,  400.,250.,+2., 120.,
     .  500.,200.,+2., 500.,  200.,200.,+2., 500.,  200.,500.,+2., 500.,
     .  500.,200.,-8., 500.,  200.,200.,-8., 500.,  200.,500.,-8., 500.,
     .  500.,200.,+8., 500.,  200.,200.,+8., 500.,  200.,500.,+8., 500.,
     .  500.,200.,-2.,1000.,  200.,200.,-2.,1000.,  200.,500.,-2.,1000.,
     .  500.,200.,+2.,1000.,  200.,200.,+2.,1000.,  200.,500.,+2.,1000.,
     .  500.,200.,-8.,1000.,  200.,200.,-8.,1000.,  200.,500.,-8.,1000.,
     .  500.,200.,+8.,1000.,  200.,200.,+8.,1000.,  200.,500.,+8.,1000./
      DATA LOU / 21 /
C
C========< Entry Point >================================================
C
C--
C  Define constants.
C--
      x2PI      = 2*ACOS(-1.)
      x4PI      = 2*x2PI
C--
C  Standard Model parameters.
C--
      S2W       = 0.23
      C2W       = 1 - S2W
      ALF       = 1/128.
      ALFS      = 0.12
      ZM        = 91.18
      WM        = 80.0
      GMZ       = 2.5
      AMZ       = ZM
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
C  sqrt(s) range.
C    NRS   = # energy points
C    RSMN  = sqrts_min [GeV] to be replaced by threshold
C            if RSMN < m_chi1 + mchi2
C    DRS   = sqrts step [GeV]
C--
      NRS  = 100
      RSMN = 400
      DRS  = 5
C--
C  Initialize ISR.
C     RSNOM = nominal sqrts [GeV] for initialization of ISR
C--
      RSNOM = 350
      CALL ISRINT(0,DBLE(RSNOM))
C--
C  Polarization.
C--
C      POL  = +0.6
C      POL  = 0.9
      POL  = 0
C      POL  = -0.9
C--
C  Couplings of electron to Z and gamma.
C--
      C0        = SQRT(x4PI*ALF)
      CW        = C0/SQRT(S2W)
      CZ        = CW/SQRT(C2W)
      QE        = -1
      T3LE      = -0.5D0
      T3RE      =  0
C--
      GAL(1)    = -C0*QE
      GAL(2)    = -C0*QE
      GZL(1)    = -CZ*(T3LE-QE*S2W)
      GZL(2)    = -CZ*(T3RE-QE*S2W)
C--
C  Prepare topdraw data.
C--
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''NEW FRAME'')')
      WRITE(LOU,'(''SET WINDOW X 2.5 12.5 Y 2.2 9.6'')')
      WRITE(LOU,'(''SET LIMITS X 200. .8E3 Y 0.001 10.'')')
      WRITE(LOU,'(''SET SCALE Y LOGARITHMIC'')')
      WRITE(LOU,'(''SET LABELS SIZE 3.0'')')
      WRITE(LOU,'(''TITLE   8. 8.8 SIZE 3.5'',
     .            '' '''' e2+3e2-3 --> C013S242+3C013S242-3'')')
      WRITE(LOU,'(''CASE                   '',
     .            '' ''''  X X X X     GUVVMVVX XGUVVMVVX X'')')
      WRITE(LOU,'(''TITLE 6.5 1.2 SIZE 4 ''''2s0O (GeV)'')')
      WRITE(LOU,'(''CASE                 ''''M UD      '')')
      WRITE(LOU,'(''TITLE 0.5 5.0 ANGLE 90. SIZE 4 '',
     .                         '' ''''S (pb)'')')
      WRITE(LOU,'(''CASE          ''''G     '')')
      WRITE(LOU,'(''SET ORDER X DUMMY Y'')')
C--
C  Loop over SUSY parameters.
C--
      DO 100 ICASE = 1, 1
C--
C  Set SUSY parameters.
C--
         AMU       = SUSY(1,ICASE,ISET)
         AM2       = SUSY(2,ICASE,ISET)
         TNB       = SUSY(3,ICASE,ISET)
         AM0       = SUSY(4,ICASE,ISET)
         AMA       = SQRT(AMU**2+AM0**2)
         BT        = ATAN(TNB)
         FM(1)     = 0.
         FM(2)     = 0.
         FM(3)     = 0.
         CALL INSUSY(AM0,AMU,AM2,TNB,AMA,ALF,ALFS,S2W,ZM,FM,
     .               SFM,SWM,SZM,GMSF,GMSW,GMSZ,
     .               GNCW,GCCZ,GNNZ,
     .               GCESNL,GCNSEL,GCNSER,
     .               GCUSDL,GCUSDR,GCDSUL,GCDSUR,
     .               GNNSNL,GNESEL,GNESER,
     .               GNUSUL,GNUSUR,GNDSDL,GNDSDR )
         AMSN      = SFM(1)
         GMX(1)    = GMSW(IXM)
         GMX(2)    = GMSW(IXP)
         GMSN      = 0
C--
C  Calculate masses and mixings.
C--
         AMX(1) = SWM(IXM)
         AMX(2) = SWM(IXP)
         QX     = -1
         IF ( IXM.NE.IXP ) QX = 0
C--
C  Coupling of chargino to Z and gamma.
C--
         GAX(1)    = -C0*QX
         GAX(2)    = -C0*QX
         GZX(1)    = GCCZ(1,IXM,IXP)
         GZX(2)    = GCCZ(2,IXM,IXP)
         GSNX(1,1) = GCESNL(1,IXM)
         GSNX(2,1) = GCESNL(2,IXM)
         GSNX(1,2) = CONJG(GCESNL(2,IXP))
         GSNX(2,2) = CONJG(GCESNL(1,IXP))
C--
C        PRINT *, ' GAX       = ', GAX
C        PRINT *, ' GZX       = ', GZX
C        PRINT *, ' GSNX(*,1) = ', GSNX(1,1), GSNX(2,1)
C        PRINT *, ' GSNX(*,2) = ', GSNX(1,2), GSNX(2,2)
C--
         SGSM = 0
         WRITE(LOU,'(''(                                    '')')
         WRITE(LOU,'(''(  **********************************'')')
         WRITE(LOU,'(''(    AMU   = '',F10.3)') SUSY(1,ICASE,ISET)
         WRITE(LOU,'(''(    AM2   = '',F10.3)') SUSY(2,ICASE,ISET)
         WRITE(LOU,'(''(    TANB  = '',F10.3)') SUSY(3,ICASE,ISET)
         WRITE(LOU,'(''(    AM0   = '',F10.3)') SUSY(4,ICASE,ISET)
         WRITE(LOU,'(''(   ---------------------------------'')')
         WRITE(LOU,'(''(    AMM   = '',F10.3)') AMX(1)
         WRITE(LOU,'(''(    AMP   = '',F10.3)') AMX(2)
         WRITE(LOU,'(''(    AMSN  = '',F10.3)') AMSN
         WRITE(LOU,'(''(  **********************************'')')
         WRITE(LOU,'(''(                                    '')')
C--
C  sqrt(s) range.
C--
         RSMN = IFIX(MAX(RSMN, AMX(1) + AMX(2) + 0.1)/5)*5 + 5
C--
         WRITE
     .      (LOU,'(''(    sqrts[GeV]     sig_0 [pb]       sig [pb]'')') 
         WRITE
     .      (LOU,'(''(   -----------------------------------------'')') 
C--
         MODE = 0
         DO 10 IRS = 0, NRS
            RS   = RSMN + IRS*DRS
            CALL SGXXEF(MODE,GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,ZM,GMZ,
     .                  AMX, GMX ,RS, POL, SG)
            CALL SGXX00(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,ZM,GMZ,
     .                  AMX, GMX ,RS, POL, SG0)
C            PRINT *, ' RS = ', RS, ' GeV  sig  = ', SG , ' pb',
C     .                                  ' sig0 = ', SG0, ' pb'
C--
            WRITE(LOU,'(3E15.5)') RS, SG0, SG
            MODE = 3
10       CONTINUE
         WRITE(LOU,'(''SET PATTERN '',A20)') PATRN(ICASE)
         WRITE(LOU,'(''JOIN '',A9)') JOIN(ICASE)
100   CONTINUE
C>>>
      PRINT *, 'That''s it, folks.'
C>>>
C--
      STOP
      END
