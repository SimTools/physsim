C* (Contents)
C*    SQRTS   : (R*8) : sqrt(s) (GeV).
C*    POLEBM  : (R*8) : electron beam polarization.
C*    ISRB    : (I*4) : ISR and BM flag.
C* (Update Record)
C*   10/11/29  K.Fujii		Original version.
C*
      COMMON /USRPRM/ SQRTS, POLEBM, ISRB
      REAL   *4       SQRTS, POLEBM
      INTEGER*4       ISRB
