CINCLUDE EI
C>>>
C*INCLUDE POTRCH
CINCLUDE POTQCD
C*INCLUDE POTCPL
C*INCLUDE POTCLM
C*INCLUDE POTMTN
C*INCLUDE POTHIG
C>>>
CINCLUDE RKUTTA3
CINCLUDE SETPRM
CINCLUDE GTGAMT
 
C
C  G_q(q) For various alpha_s
C
C
      IMPLICIT REAL*8 ( A-H, O-Z )
      COMPLEX*16      GQ
      CHARACTER*8 JOIN(0:3)
      DATA JOIN / 'DASH', 'SOLID', 'DOTDASH', 'DOT' /
C--
C  Initialize parameters.
C--
      PI   = ACOS(-1.D0)
      ALP  = 1./128.
      SN2W = 0.23
      AMSZ = 91.17
      GAMZ = 2.5
      AMSW = 80.0
      AMSB = 5.0
      AMST = 150
      VFF  = 1.0
      AMSH = 1.D8
      BETH = 0
C--
C  Alpha_s range.
C--
      NA     = 2
      ALFSMN = 0.11
      ALFSMX = 0.13
      DA     = (ALFSMX-ALFSMN)/NA
C>>>
      NA     = 0
      ALFSMN = 0.12
C     NA     = 0
C     ALFSMN = 0.16
C--
C  Energy range.
C--
C     NE   = 2
C     EMIN =-2
C     EMAX = 2
C     DE   = (EMAX-EMIN)/NE
C>>>
      NE   = 3
      EMIN =-6
      EMAX = 3
      DE   = (EMAX-EMIN)/NE
C>>>
      NE   = 0
      EMIN = 0
C--
C  Momentum range.
C--
      NQ   = 50
      QMAX = 50
C     QMAX  = SQRT( (3*AMST+EMIN-AMSW)*(3*AMST+EMIN+AMSW)
C    .           *(AMST+EMIN-AMSW)*(AMST+EMIN+AMSW) )/2/(2*AMST+EMIN)
      DQ   = QMAX/NQ
C--
C  Start drawing G_q(q).
C--
      DO 10000 IA = 0, NA
         ALPS = ALFSMN + IA*DA
         WRITE(20,*) '( ALFS = ', ALPS
         MODE = 0
         DO 1000 IE = 0, NE
            E = EMIN + IE*DE
            RSH = 2*AMST + E
            WRITE(20,*) '( E = ', E, ' GeV'
            WRITE(20,*) 'SET ORDER X Y'
            DO 100 IQ = 0, NQ
               Q    = IQ*DQ
               CALL GREENQ(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                       AMST, VFF, AMSH,BETH,   RSH,Q,GQ)
C>>>
C              WRITE(20,*) Q, ABS(GQ)**2
               WRITE(20,*) Q, 4*PI*ABS(GQ)**2*Q**2
C>>>
               MODE = 3
100         CONTINUE
            WRITE(20,*) 'JOIN ', JOIN(IA)
            MODE = 2
1000     CONTINUE
10000 CONTINUE
C--
C  That's it.
C--
      END
