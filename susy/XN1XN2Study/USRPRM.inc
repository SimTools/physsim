C* (Contents)
C*    SQRTS   : (R*8) : sqrt(s).
C*    POLEBM  : (R*8) : electron beam polarization.
C*    SGMEBM  : (R*8) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C* (Update Record)
C*   95/03/21  K.Fujii		Original version.
C*

      COMMON /USRPRM/ SQRTS, POLEBM, SGMEBM, GAMSZ2, ISRB
      REAL   *4       SQRTS, POLEBM, SGMEBM, GAMSZ2
      INTEGER*4       ISRB
