C* (Contents)
C*    SQRTS   : (R*4) : sqrt(s).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IMD1LO  : (I*4) : 1st W- decay mode ID.
C*    IMD1HI  : (I*4) : lst W- decay mode ID.
C*    IMD2LO  : (I*4) : 1st W+ decay mode ID.
C*    IMD2HI  : (I*4) : lst W+ decay mode ID.
C* (Update Record)
C*   95/03/21  K.Fujii		Original version.
C*   09/09/14  T.Tanabe         Added alternative CS switch.
C*
      COMMON /USRPRM/ SQRTS, POLEBM, SGMEBM, ISRB,
     .                IMD1LO, IMD1HI, IMD2LO, IMD2HI,
     .                POLPBM, ALTCS
#ifdef WITH_DBD_STANDARD
     .              ,  NDCYSL, ISR_LLA_ORDER, 
     .                 NS_REMNANTS, NS_BEAMS
      INTEGER*4        NDCYSL, ISR_LLA_Order
      INTEGER*4        NS_REMNANTS, NS_BEAMS
#endif

      REAL   *4       SQRTS, POLEBM, SGMEBM, POLPBM
      INTEGER*4       ISRB,  
     .                IMD1LO, IMD1HI, IMD2LO, IMD2HI,
     .                ALTCS
