C* (Contents)
C*    SQRTS   : (R*4) : sqrt(s).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IGFR    : (I*4) : gneration number.
C*    ITFR    : (I*4) : isospin, (1,2)=(up,down).
C*    LQFR    : (I*4) : lepton-quark switch, (1,2)=(l,q).
C*    PECT    : (R*4) : e- energy cut (GeV).
C*    CECL    : (R*4) : cos(theta_e-)_min
C*    CECU    : (R*4) : cos(theta_e-)_max
C*    PPCT    : (R*4) : e+ energy cut (GeV).
C*    CPCL    : (R*4) : cos(theta_e+)_min
C*    CPCU    : (R*4) : cos(theta_e+)_max
C*    PFCT    : (R*4) : fermion energy cut (GeV).
C*    CFCT    : (R*4) : cos(theta_f) cut.
C*    WMNF    : (R*4) : minimum m_ff (GeV).
C* (Update Record)
C*   96/03/07  K.Fujii		Original version.
C*

      COMMON /USRPRM/ SQRTS, POLEBM, SGMEBM, ISRB, 
     .                IGFR, ITFR, LQFR, PECT, CECL, CECU, 
     .                PPCT, CPCL, CPCU, PFCT, CFCT, WMNF
      REAL   *4       SQRTS, POLEBM, SGMEBM, PECT, CECL, CECU,
     .                PPCT, CPCL, CPCU, PFCT, CFCT, WMNF
      INTEGER*4       ISRB, IGFR, ITFR, LQFR
