C* (Contents)
C*    DELTRS  : (R*4) : sqrt(s) - 2 m_t (GeV).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IMD1LO  : (I*4) : 1st W- decay mode ID.
C*    IMD1HI  : (I*4) : lst W- decay mode ID.
C*    IMD2LO  : (I*4) : 1st W+ decay mode ID.
C*    IMD2HI  : (I*4) : lst W+ decay mode ID.
C*    BETH    : (R*4) : top Yukawa normalized by SM.
C*    VKMT    : (R*4) : V_KM for top normalized by SM.
C*    NQCD    : (I*4) : (0,1)=(w/o NRQCD, w/ NRQCD).
C* (Update Record)
C*   95/03/16  K.Fujii		Original version.
C*   99/09/12  K.Fujii		Added W's decay mode selectors.
C*
      COMMON /USRPRM/ DELTRS, POLEBM, SGMEBM, ISRB,
     .                IMD1LO, IMD1HI, IMD2LO, IMD2HI,
     .                VKMT, BETH, NQCD
      REAL   *4       DELTRS, POLEBM, SGMEBM, VKMT, BETH
      INTEGER*4       ISRB, IMD1LO, IMD1HI, IMD2LO, IMD2HI, NQCD
