C* (Contents)
C*    SQRTS   : (R*8) : sqrt(s) (GeV).
C*    POLEBM  : (R*8) : electron beam polarization.
C*    ISRB    : (I*4) : ISR and BM flag.
C*    IMD1LO  : (I*4) : 1st W- decay mode ID.
C*    IMD1HI  : (I*4) : lst W- decay mode ID.
C*    IMD2LO  : (I*4) : 1st W+ decay mode ID.
C*    IMD2HI  : (I*4) : lst W+ decay mode ID.
C* (Update Record)
C*   10/11/29  K.Fujii		Original version.
C*
      COMMON /USRPRM/ SQRTS, POLEBM, ISRB,
     .                IMD1LO, IMD1HI, IMD2LO, IMD2HI
      REAL   *4       SQRTS, POLEBM
      INTEGER*4       ISRB, IMD1LO, IMD1HI, IMD2LO, IMD2HI
