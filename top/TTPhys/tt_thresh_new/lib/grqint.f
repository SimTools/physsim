 
      SUBROUTINE  GRQINT(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,AMST,VTB2,E)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
      INTEGER*4       MODE
      REAL*8          ALPS, ALP, SN2W, AMSZ, AMSW, AMSB, AMST, VTB2, E
C--
C  Parameters.
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /RKPARM/ DRS, DRP, EPSS, EPSP
      REAL   *8       DRS, DRP, EPSS, EPSP
C--
      COMMON /QCDPAR/ ALFST, ALFSP, RQCD
      REAL*8          ALFST, ALFSP, RQCD
C--
C  Other variables.
C--
      REAL   *8       GAMAT0, GF, RADCOR, H, KAP2
C--
      REAL   *8       GMTWMB
C--
C  Block data.
C--
      EXTERNAL        GRQBLK
      SAVE
C
C===========< Entry Point >=============================================
C
C--
C  Set up constants.
C     MODE = 0 : 1-st call.
C          = 1 : if alpha_s changed.
C          = 2 : some parameters but alpha_s changed.
C          = 3 : nothing but E changed.
C--
      IF ( MODE.EQ.0 ) THEN
         PI   = ACOS(-1.D0)
         IMAG = (0.D0,1.D0)
         MB   = AMSB
         MZ   = AMSZ
         MW   = AMSW
         CF   = 4.D0/3
         ALF  = ALP
         S2W  = SN2W
         GF   = SQRT(2.D0)*4*PI*ALP/SN2W/8/AMSW**2
      ENDIF
C--
C  Initialize variables.
C--
      IF ( MODE.LE.2 ) THEN
         ALFS = ALPS
         MT   = AMST
C--
C  Set y_cut.
C--
         YCUT = (1-MW/MT)**2 * .999999D0
C>>>
         PRINT *, 'GRQINT sets y_cut at y_max.'
C>>>
C--
C  Calculate Lambda_cutoff.
C--
         ALAMB = SQRT( (9*MT**2-MW**2)*(MT**2-MW**2) )/4/MT
C--
C  Calculate Gamma_t0.
C--
         RWT    = MW**2/MT**2
         KAPPA  = ( 1 - 2*RWT )/( 1 + 2*RWT )
         GAMAT0 = GF/SQRT(2.D0) * MT**3/8/PI * VTB2
     .               * ( 1 + 2*RWT )*( 1 - RWT )**2
C--
C  Initialize QCD potential.
C--
         IF ( MODE.LE.1 ) THEN
            CALL SETVPR(ALFS,MT,MB,MZ,GAMAT0,ALFST,ALFSP,RQCD)
            PRINT *, 'alpha_s(mZ)      = ', ALFS
            PRINT *, 'alpha_s(mt)      = ', ALFST
            PRINT *, 'alpha_s(mt*alfs) = ', ALFSP
         END IF
C--
C  QCD correction to Gamma_t0.
C--
         RADCOR = 1 - CF*ALFST/2/PI * H(RWT)
         GAMMAT = GAMAT0 * RADCOR
C>>>
         WRITE(*,*)  'Width without RADCOR = ', GAMAT0
         WRITE(*,*)  'Width with    RADCOR = ', GAMMAT
C>>>
         WRITE(*,*)  'Ratio of the width with b mass',
     .                GMTWMB( GF, MT, MW, MB )/GAMAT0
C        WRITE(*,*)  'check',
C    .                GMTWMB( GF, MT, MW, 0.D0 )/GAMAT0
C>>>
C   The previous treatment included
C      PS suppression + time dilatation
C   but Coulombic enhancement.
C
C1       TLETA1 = 3*( 1 + RWT + 2*RWT**2 ) /( 1 + RWT - 2*RWT**2 )
C1       TLETA2 = 1/6.D0*( 13 + 13*RWT + 46*RWT**2 )
C1   .                        /( 1 + RWT - 2*RWT**2 )
C
C 94/08/01 New Running Width (TUW-94-06).
C   In the case of an exactly Coulombic potential,
C      PS suppression + Coulombic enhancement = 0
C   which leaves only time dilatation effect.
C
         PRINT *, 'Running width w/ time dilatation only.'
         TLETA1 = 0
         TLETA2 = 1/2.D0
C>>>
         KAP2   = GAMMAT*TLETA2/MT
         EFFMT  = MT/( 1 + IMAG*KAP2 )
      END IF
C--
C  If the tabulation of V(r) is not necessary.
C--
      IF ( MODE.LE.2 ) THEN
C--
C  Loop controlling variables.
C--
C>>>
         DRS   = 1.D-4
         DRP   = 1.D-4
C        DRS   = 1.D-3
C        DRP   = 1.D-3
C>>>
         EPSS  = 1.D-13
         EPSP  = 1.D-13
C--
C  Set up parameters related to the QCD potential.
C--
         CALL GTQCDP( RQCD )
      END IF
C--
C  Normalization of T0, T1, T2.
C--
      CALL GTCEFF(E)
C--
C  That's it.
C--
      RETURN
      END
