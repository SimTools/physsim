      block data blkdat
C--
C  THIS FILE IS AUTOMATICALLY GENERAGED FROM
C     include/SMPTAB.inc
C  DO NOT EDIT. EDIT C= COMMENTS IN THESE INCLUDE FILES
C--
C*
C*==================
C* Common /SMPTAB/
C*==================
C*
C*    This common contains parton properties such as masses and
C*    widths of elementary fermions and gauge bosons.
C*    The constants in this common are initialized by CONINT.
C*
C*
      COMMON /SMPTAB/ AMF(3,2,2), GMF(3,2,2), AMW, AMZ, AMH,
     .                GMZTOT, GMWTOT, GMHTOT, MDVDK(3,12,2),
     .                BRVDK(0:12,2), VKM(3,3,2)
      INTEGER*4       MDVDK
      REAL   *8       AMF, GMF, AMW, AMZ, AMH,
     .                GMZTOT, GMWTOT, GMHTOT,
     .                BRVDK, VKM
C--
C  Fermion mass (IG,IT,LQ).
C--
      DATA AMF   /  1.D-10      ,   1.D-10    ,   1.D-10  ,
     .              0.5110034D-3, 105.65946D-3,   1.7770D0,
     .              0.04D0      ,   1.5D0     , 175.0D0   ,
     .              0.04D0      ,   0.1D0     ,   4.7D0   /
C--
C  Fermion width (IG,IT,LQ).
C--
      DATA GMF   /  1.D-10      ,   1.D-10    ,   1.0D-10 ,
     .              1.D-10      ,   1.D-10    ,   1.0D-10 ,
     .              1.D-10      ,   1.D-10    ,   0.2D+01 ,
     .              1.D-10      ,   1.D-10    ,   1.0D-10 /
C--
C  Gauge boson and Higgs mass and width.
C--
      DATA AMZ   /  91.1D0  /   GMZTOT /  2.3541D0 /
      DATA AMW   /  80.0D0  /   GMWTOT /  1.9127D0 /
      DATA AMH   /  90.0D0  /   GMHTOT /  1.0000D0 /
C--
C  KM matrix.
C--                   nu_e           nu_mu         nu_tau
      DATA VKM   /  1.000D0     ,   0.000D0   ,   0.000D0,
     .              0.000D0     ,   1.000D0   ,   0.000D0,
     .              0.000D0     ,   0.000D0   ,   1.000D0,
C                      u               c             t
     .              0.975D0     ,   0.222D0   ,   0.010D0,
     .              0.222D0     ,   0.974D0   ,   0.043D0,
     .              0.010D0     ,   0.043D0   ,   0.999D0/
C--
C  W decay mode (IU,ID,LQ)
C--
      DATA MDVDK  / 1,1,1, 2,2,1, 3,3,1,  1,1,2, 2,1,2, 3,1,2,
     .              1,2,2, 2,2,2, 3,2,2,  1,3,2, 2,3,2, 3,3,2,
C--
C  Z decay mode (IG,IT,LQ)
C--
     .              1,1,1, 2,1,1, 3,1,1,  1,2,1, 2,2,1, 3,2,1,
     .              1,1,2, 2,1,2, 3,1,2,  1,2,2, 2,2,2, 3,2,2/
C--
C  Vector boson decay branching fractions.
C--
      DATA BRVDK  / 26*0.D0 /
C--
C  THIS FILE IS AUTOMATICALLY GENERAGED FROM
C      include/BEMCNS.inc
C  DO NOT EDIT. EDIT C= COMMENTS IN THESE INCLUDE FILES
C--
      COMMON /BEMCNS/ SGEBM
      REAL*8          SGEBM
      DATA SGEBM / 0.007 /
      end
