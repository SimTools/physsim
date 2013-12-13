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
      COMPLEX*16      GQ, GQT
      CHARACTER*8 JOIN(0:3)
      DATA JOIN / 'DASH', 'SOLID', 'DOTDASH', 'DOT' /
C--
C  Initialize parameters.
C--
      PI   = ACOS(-1.D0)
      ALP  = 1./128.
      SN2W = 0.23
      AMSZ = 91.18
      GAMZ = 2.5
      AMSW = 80.0
      AMSB = 5.0
      AMST = 175
      VFF  = 1.0
      AMSH = 1.D8
      BETH = 0
#if 1
      xGV2FB = 3.8937966E11
      S2W   = SN2W
      SCW   = SQRT(S2W*(1-S2W))
      QE    = -1
      VE    = ( -1.D0/2 - 2*QE*S2W )/2/SCW
      AE    = ( -1.D0/2            )/2/SCW
      QT    = 2.D0/3
      VT    = ( +1.D0/2 - 2*QT*S2W )/2/SCW
#endif
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
#if 1
      NE   = 70
      EMIN =-20
      EMAX = 50
      DE   = (EMAX-EMIN)/NE
#else
      NE   = 0
      EMIN = 0
      EMAX = 0
      DE   = 0
#endif
C--
C  Start drawing G_q(q).
C--
      DO 10000 IA = 0, NA
         ALPS = ALFSMN + IA*DA
         WRITE(20,*) '( ALFS = ', ALPS
         MODE = 0
#if 1
         WRITE(21,*) '( alps = ', alps
         WRITE(21,*) 'SET ORDER X Y'
#endif
         DO 1000 IE = 0, NE
            E = EMIN + IE*DE
            RSH = 2*AMST + E
C--
C  Momentum range.
C--
#if 0
            NQ   = 50
            QMAX = 50
#else
            NQ   = 500
            QMAX = SQRT( (3*AMST+E-AMSW)*(3*AMST+E+AMSW)
     .           *(AMST+E-AMSW)*(AMST+E+AMSW) )/2/(2*AMST+E)
#endif
            DQ   = QMAX/NQ
#if 1
            sum = 0
#endif
            WRITE(20,*) '( E = ', E, ' GeV'
            WRITE(20,*) 'SET ORDER X Y'
            DO 100 IQ = 0, NQ
               Q    = IQ*DQ
#if 0
               CALL GREENQ(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                       AMST, VFF, AMSH,BETH,   RSH,Q,GQ)
#else
               call gqtabl(alps,amst,amsh,beth,rsh,q,gqt)
#if 0
               CALL GREENQ(MODE,ALPS,ALP,SN2W,AMSZ,GAMZ,AMSW,AMSB,
     .                       AMST, VFF, AMSH,BETH,   RSH,Q,GQ)
               print *, ' rsh = ', rsh, ' q = ', q
               print *, ' gqtab = ', gqt
               print *, ' gq    = ', gq
#endif
               gq = gqt
#endif
               IF ( IQ.EQ.0 .OR. IQ.EQ.NQ ) THEN
                  WT = 1
               ELSE IF ( MOD(IQ,2).EQ.0 ) THEN
                  WT = 2
               ELSE
                  WT = 4
               ENDIF
               add = 4*PI*ABS(GQ)**2*Q**2
               WRITE(20,*) Q, add
               call gmthet(q,gm)
               sum = sum +  add*gm*dq*wt/3
               MODE = 3
100         CONTINUE
            WRITE(20,*) 'JOIN ', JOIN(IA)
            MODE = 2
#if 1
            S      = RSH**2
            DZ     = S/((RSH-AMSZ)*(RSH+AMSZ))
            SG     = (6*(ALP)**2/PI)
     .         *( QT**2 - 2*QT*VE*VT*DZ + (VE**2+AE**2)*VT**2*DZ**2)
     .         *( 1 - 8*ALPS/3/PI )**2/S**2
     .         *xGV2FB*sum
            print *, rsh, sg
            WRITE(21,*) rsh, sg
#endif
1000     CONTINUE
#if 1
         WRITE(21,*) 'JOIN ', JOIN(IA)
#endif
10000 CONTINUE
C--
C  That's it.
C--
      END
