CC**********************************************************************
C*
C*===================================------===
C* Subroutine UHQ2BW(AM,GM,QMN,QMX,Z,Q2,WAT)
C*===================================------===
C*
C* (Purpose)
C*    Reparametrizes a Breit-Wigner resonance for less singular
C*    behavior in numerical integration.
C* (Inputs)
C*       AM     : (R*16) : mass.
C*       GM     : (R*16) : width.
C*       QMN    : (R*16) : Q_min.
C*       QMX    : (R*16) : Q_max.
C*       Z      : (R*16) : integration variable in (0,1).
C* (Output)
C*       Q2     : (R*16) : Q^2.
C*       WAT    : (R*16) : Jacobian weight.
C* (Relation)
C*    Invokes no external routines.
C* (Update Record)
C*   95/04/07  K.Fujii	Original version.
C*
CC**********************************************************************

      SUBROUTINE UHQ2BW(AM,GM,QMN,QMX,Z,Q2,WAT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      REAL   *8  AM, GM, QMN, QMX, Z, Q2, WAT
C
C========< Entry Point >================================================
C
C--
C  Define some variables.
C--
      AMGM  = AM*GM
      AM2   = AM*AM
      AMGM2 = AMGM*AMGM
C--
C  Calculate Q2.
C--
      THMN = ATAN((QMN-AM)*(QMN+AM)/AMGM)
      THMX = ATAN((QMX-AM)*(QMX+AM)/AMGM)
      TH   = THMN + (THMX-THMN)*Z
      Q2   = AMGM*TAN(TH) + AM2
C--
C  Then calculate Jacobian.
C--
      WAT  = (THMX-THMN)*( (Q2-AM2)**2 + AMGM2 )/AMGM
#if 0
C
C-- Narrow width approximation.      
C
      if ( AM .GE. QMN .AND. AM.LE. QMX ) THEN
         xPI  = ACOS(-1.D0)
         Q2   = AM2
         WAT  = xPI*AMGM
      ENDIF
#endif
C--
C  That's it.
C--
      RETURN
      END
