CC**********************************************************************
C*
C*============================-----------===
C* Subroutine GETXPH(AM,AMR,Z,IPV,XPH,WAT)
C*============================-----------===
C*
C* (Purpose)
C*    Generate phase space variables for a given set of random numbers
C*    of this iteration.
C* (Inputs)
C*       AM(i)    : (R*4) : mass  of i.
C*       AMR(1,r) : (R*4) : mass of resonance r        
C*          (2,r) : (R*4) : width of resonance r           
C*          (3,r) : (R*4) : m_(ij) + m_k  for r <--> (ij)
C*       Z(1-7)   : (R*8) : integration varialble in (0,1).
C* (Outputs)
C*       IPV(1)   : (I*4) : a+
C*          (2)   : (I*4) : b+
C*          (3)   : (I*4) : c+
C*       XPH(1,1) : (R*4) : invariant mass squared for X+.
C*          (2,1) : (R*4) : cos(theta_a+) in X+ rest frame.
C*          (3,1) : (R*4) : phi_a+ in X+ rest frame.
C*       XPH(1,2) : (R*4) : invariant mass squared for a+b+.
C*          (2,2) : (R*4) : invariant mass squared for a+c+.
C*          (3,2) : (R*4) : phi_b+ in X+ rest frame.
C*       WAT      : (R*8) : Jacobian weight.
C* (Relation)
C*    Invokes UHQ2BW and UHQIJ3 in ttz_lib.
C* (Update Record)
C*   95/04/10  K.Fujii	Original version.
C*
CC**********************************************************************

      SUBROUTINE GETXPH(AM,AMR,Z,IPV,XPH,WAT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      INTEGER*4  IPV(*)
      REAL   *4  AM(*), AMR(3,3), XPH(3,*)
      REAL   *8  Z(*), WAT
C--
C  Local variables.
C--
      INTEGER*4  JPV(3,3,3)
      REAL   *8  WT
      DATA EPS   / 5.D0 /
      DATA NCALL /    0 /
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
C--
         x2PI  = 2*ACOS(-1.D0)
C--
         JPV(1,1,2) = 2
         JPV(2,1,2) = 3
         JPV(3,1,2) = 4
         JPV(1,1,3) = 3
         JPV(2,1,3) = 2
         JPV(3,1,3) = 4
         JPV(1,2,1) = 2
         JPV(2,2,1) = 4
         JPV(3,2,1) = 3
         JPV(1,2,3) = 4
         JPV(2,2,3) = 2
         JPV(3,2,3) = 3
         JPV(1,3,1) = 3
         JPV(2,3,1) = 4
         JPV(3,3,1) = 2
         JPV(1,3,2) = 4
         JPV(2,3,2) = 3
         JPV(3,3,2) = 2
      ENDIF
C--
C  Set independent variables.
C     Z( 1) : m(+ab)**2
C      ( 2) : m(+ac)**2
C      ( 3) : cos(theta_a)     in X+ rest frame
C      ( 4) : phi_a            in X+ rest frame
C      ( 5) : phi_b            in X+ rest frame
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Set invariant masses squared for daughter particles from X+.
C--
      AMX = AM(1)
      CALL UHQIJ3(AM(1),JPV,AMR(1,1),EPS,AMX,Z(1),IPV(1),Q42,Q52,WT)
      WAT = WAT*WT
C>>>
CCDBG       print *, ' IPV(*) = ', (IPV(K),K=1,3)
C>>>       
C--
C  Add angular variables.
C--
      XPH(1,1) = AMX*AMX
      XPH(2,1) = -1 + 2*Z(3)
      XPH(3,1) = x2PI*Z(4)
      XPH(1,2) = Q42
      XPH(2,2) = Q52
      XPH(3,2) = x2PI*Z(5)
C--
C  Calculate Jacobian.
C--
      WAT = WAT*2*(x2PI)**2
C--
C  That's it.
C--
      RETURN
      END
