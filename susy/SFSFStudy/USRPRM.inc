C* (Contents)
C*    SQRTS   : (R*4) : sqrt(s).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C*    IGNSF   : (R*4) : generation number.
C*    IHNDM   : (R*4) : (1,2) for (L,R) for sf-.
C*    IHNDP   : (R*4) : (1,2) for (L,R) for sf+.
C* (Update Record)
C*   95/03/21  K.Fujii		Original version.
C*

      COMMON /USRPRM/ SQRTS, POLEBM, SGMEBM, IGNSF, IHNDM, IHNDP,
     .                IDOTU, HTUM, HTUP, ISRB
      REAL   *4       SQRTS, POLEBM, SGMEBM, HTUM, HTUP
      INTEGER*4       IGNSF, IHNDM, IHNDP, IDOTU, ISRB
