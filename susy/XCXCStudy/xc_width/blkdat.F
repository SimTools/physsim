      block data blkdat
C--
C  THIS FILE IS AUTOMATICALLY GENERAGED FROM
C     include/SMPTAB.inc
C     include/USMPRM.inc
C     include/BEMCNS.inc
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
      REAL   *4       AMF, GMF, AMW, AMZ, AMH,
     .                GMZTOT, GMWTOT, GMHTOT,
     .                BRVDK, VKM
C--
C  Fermion mass (IG,IT,LQ).
C--
      DATA AMF   /  1.E-10      ,   1.E-10    ,   1.E-10 ,
     .              0.5110034E-3, 105.65946E-3,   1.732E0,
     .              0.04E0      ,   1.5E0     , 135.0E0  ,
     .              0.04E0      ,   0.1E0     ,   4.7E0  /
C--
C  Fermion width (IG,IT,LQ).
C--
      DATA GMF   /  1.E-10      ,   1.E-10    ,   1.E-10 ,
     .              1.E-10      ,   1.E-10    ,   1.E-10 ,
     .              1.E-10      ,   1.E-10    ,   1.E-10 ,
     .              1.E-10      ,   1.E-10    ,   1.E-10 /
C--
C  Gauge boson and Higgs mass and width.
C--
      DATA AMZ   /  91.1E0  /   GMZTOT /  2.3541E0 /
      DATA AMW   /  80.0E0  /   GMWTOT /  1.9127E0 /
      DATA AMH   /  90.0E0  /   GMHTOT /  1.0000E0 /
C--
C  KM matrix.
C--                   nu_e           nu_mu         nu_tau
      DATA VKM   /  1.000E0     ,   0.000E0   ,   0.000E0,
     .              0.000E0     ,   1.000E0   ,   0.000E0,
     .              0.000E0     ,   0.000E0   ,   1.000E0,
C                      u               c             t
     .              0.975E0     ,   0.222E0   ,   0.010E0,
     .              0.222E0     ,   0.974E0   ,   0.043E0,
     .              0.010E0     ,   0.043E0   ,   0.999E0/
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
      DATA BRVDK  / 26*0. /
C* (Contents)
C*    ALFI  : (R*4) : 1/alpha(m_Z).
C*    ALFS  : (R*4) : alpha_s(m_Z).
C*    AMSW  : (R*4) : m_W (GeV).
C*    AMSZ  : (R*4) : m_Z (GeV).
C*    AMSH  : (R*4) : m_H (GeV).
C*    AMST  : (R*4) : m_t (GeV).
C*
C* (Update Record)
C*   95/03/16  K.Fujii		Original version.
C*
      COMMON /USMPRM/  ALFI, ALFS, AMSW, AMSZ, AMSH, AMST
      REAL*4           ALFI, ALFS, AMSW, AMSZ, AMSH, AMST
      DATA    ALFI  /  128.0 /
      DATA    ALFS  /  0.120 /
      DATA    AMSW  /  80.00 /
      DATA    AMSZ  /  91.17 /
      DATA    AMSH  /  300.0 /
      DATA    AMST  /  170.0 /
      COMMON /BEMCNS/ SGEBM
      REAL*4          SGEBM
      DATA SGEBM / 0.007 /
      end
