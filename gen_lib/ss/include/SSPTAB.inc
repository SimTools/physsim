C--
C  1st and 2nd Generations
C     SFM(I): (SNL,SEL,SER,SUL,SUR,SDL,SDR)
C     SHM   : (AMHL,AMHH,AMA,AMHC)
C     ALPH  : HIGGS MIXING ANGLE
C--
      COMMON /SSPTAB/ SWM, SZM, SFM, SHM, GMSW, GMSZ, GMSF, GMSH, ALPH
      REAL*4          SWM(2), SZM(4), SFM(7), SHM(4), 
     .                GMSW(2), GMSZ(4), GMSF(7), GMSH(4)
      REAL*4 ALPH
C--
C  3rd Generation
C     SFM3M(I): masses for (SNL,SEL,SER,SUL,SUR,SDL,SDR)
C     SPHI(IG): mixing angles for (STAU,STOP,SBOTTOM) 
C--
      COMMON /SSPTAB2/ SFM3M, SPHI, GMSF3
      REAL*4           SFM3M(7), SPHI(3), GMSF3(7)
