C* (Contents)
C*    SQRTS   : (R*4) : sqrt(s).
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IBTP    : (I*4) : beam parameter flag.
C* (Update Record)
C*   95/06/23  K.Fujii		Original version.
C*
      COMMON /USRPRM/ SQRTS, SGMEBM, ISRB, IBTP
      REAL   *4       SQRTS, SGMEBM
      INTEGER*4       ISRB, IBTP
