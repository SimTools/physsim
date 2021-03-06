C* (Contents)
C*    SQRTS   : (R*4) : sqrt(s).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IMD1LO  : (I*4) : 1st W+ decay mode ID.
C*    IMD1HI  : (I*4) : lst W+ decay mode ID.
C*    IMD2LO  : (I*4) : 1st Z  decay mode ID.
C*    IMD2HI  : (I*4) : lst Z  decay mode ID.
C* (Update Record)
C*   95/03/21  K.Fujii		Original version.
C*

      COMMON /USRPRM/ SQRTS, POLEBM, SGMEBM, ISRB,
     .                IMD1LO, IMD1HI, IMD2LO, IMD2HI
      REAL   *4       SQRTS, POLEBM, SGMEBM
      INTEGER*4       ISRB,  IMD1LO, IMD1HI, IMD2LO, IMD2HI
