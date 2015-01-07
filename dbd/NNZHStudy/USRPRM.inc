C* (Contents)
C*    SQRTS   : (R*8) : sqrt(s).
C*    POLEBM  : (R*8) : electron beam polarization.
C*    SGMEBM  : (R*8) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C* (Update Record)
C*   95/03/21  K.Fujii		Original version.
C*

      COMMON /USRPRM/ SQRTS, POLEBM, SGMEBM, ISRB, IMDMIN, IMDMAX
#ifdef WITH_DBD_STANDARD
     .              ,  ISR_LLA_ORDER, 
     .                 NS_REMNANTS, NS_BEAMS
      INTEGER*4        ISR_LLA_Order
      INTEGER*4        NS_REMNANTS, NS_BEAMS
#endif
      REAL   *4       SQRTS, POLEBM, SGMEBM
      INTEGER*4       ISRB,  IMDMIN, IMDMAX
