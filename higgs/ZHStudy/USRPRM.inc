C* (Contents)
C*    SQRTS   : (R*4) : sqrt(s).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IMD1LO  : (I*4) : minimum decay mode ID for Z0.
C*    IMD1HI  : (I*4) : maximum decay mode ID for Z0.
C*    IMD2LO  : (I*4) : minimum decay mode ID for H0 (ffbar only).
C*    IMD2HI  : (I*4) : maximum decay mode ID for H0 (ffbar only).
C* (Update Record)
C*   95/03/21  K.Fujii		Original version.
C*   99/05/21  K.Fujii		Added decay mode control.
C*

      COMMON /USRPRM/ SQRTS,  POLEBM, SGMEBM, ISRB, 
     .                IMD1LO, IMD1HI, IMD2LO, IMD2HI
      REAL   *4       SQRTS, POLEBM, SGMEBM
      INTEGER*4       ISRB,  IMD1LO, IMD1HI, IMD2LO, IMD2HI
