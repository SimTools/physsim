C* (Contents)
C*    DELTRS  : (R*4) : sqrt(s) - 2 m_t (GeV).
C*    POLEBM  : (R*4) : electron beam polarization.
C*    SGMEBM  : (R*4) : beam energy spread (fraction).
C* (Update Record)
C*   95/03/16  K.Fujii		Original version.
C*
      COMMON /USRPRM/ DELTRS, POLEBM, SGMEBM
      REAL*4          DELTRS, POLEBM, SGMEBM
C
C=    DATA  DELTRS  / -0.595 /
C=    DATA  POLEBM  /  0.0   /
C=    DATA  SGMEBM  /  0.005 /
