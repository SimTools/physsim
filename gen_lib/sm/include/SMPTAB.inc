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
      REAL   *4       AMF, GMF, AMW, AMZ, AMH,
     .                GMZTOT, GMWTOT, GMHTOT,
     .                BRVDK, VKM
C--
C  Fermion mass (IG,IT,LQ).
C--
C=    DATA AMF   /  1.E-10      ,   1.E-10    ,   1.E-10 ,
C=   .              0.5110034E-3, 105.65946E-3,   1.732E0,
C=   .              0.04E0      ,   1.5E0     , 135.0E0  ,
C=   .              0.04E0      ,   0.1E0     ,   4.7E0  /
C--
C  Fermion width (IG,IT,LQ).
C--
C=    DATA GMF   /  1.E-10      ,   1.E-10    ,   1.E-10 ,
C=   .              1.E-10      ,   1.E-10    ,   1.E-10 ,
C=   .              1.E-10      ,   1.E-10    ,   1.E-10 ,
C=   .              1.E-10      ,   1.E-10    ,   1.E-10 /
C--
C  Gauge boson and Higgs mass and width.
C--
C=    DATA AMZ   /  91.1E0  /   GMZTOT /  2.3541E0 /
C=    DATA AMW   /  80.0E0  /   GMWTOT /  1.9127E0 /
C=    DATA AMH   /  90.0E0  /   GMHTOT /  1.0000E0 /
C--
C  KM matrix.
C--                   nu_e           nu_mu         nu_tau
C=    DATA VKM   /  1.000E0     ,   0.000E0   ,   0.000E0,
C=   .              0.000E0     ,   1.000E0   ,   0.000E0,
C=   .              0.000E0     ,   0.000E0   ,   1.000E0,
C                      u               c             t
C=   .              0.975E0     ,   0.222E0   ,   0.010E0,
C=   .              0.222E0     ,   0.974E0   ,   0.043E0,
C=   .              0.010E0     ,   0.043E0   ,   0.999E0/
C--
C  W decay mode (IU,ID,LQ)
C--
C=    DATA MDVDK  / 1,1,1, 2,2,1, 3,3,1,  1,1,2, 2,1,2, 3,1,2,
C=   .              1,2,2, 2,2,2, 3,2,2,  1,3,2, 2,3,2, 3,3,2,
C--
C  Z decay mode (IG,IT,LQ)
C--
C=   .              1,1,1, 2,1,1, 3,1,1,  1,2,1, 2,2,1, 3,2,1,
C=   .              1,1,2, 2,1,2, 3,1,2,  1,2,2, 2,2,2, 3,2,2/
C--
C  Vector boson decay branching fractions.
C--
C=    DATA BRVDK  / 26*0. /
