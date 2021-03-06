CC**********************************************************************
C*
C*=======================================-----------===
C* Subroutine GETXPH(RS,AM,AMX,GMX,AMR,Z,IPV,XPH,WAT)
C*=======================================-----------===
C*
C* (Purpose)
C*    Generate phase space variables for a given set of random numbers
C*    of this iteration.
C* (Inputs)
C*       RS         : (R*4) : sqrt(s).
C*       AM(i)      : (R*4) : mass  of i.
C*       AMX(1)     : (R*4) : mass  of X-.
C*          (2)     : (R*4) : mass  of X+.
C*       GMX(1)     : (R*4) : width of X-.
C*          (2)     : (R*4) : width of X+.
C*       AMR(1,r,1) : (R*4) : mass of resonance r        
C*          (2,r,1) : (R*4) : width of resonance r           for X+.
C*          (3,r,1) : (R*4) : m_(ij) + m_k  for r <--> (ij)
C*       AMR(1,r,2) : (R*4) : mass of resonance r        
C*          (2,r,2) : (R*4) : width of resonance r           for X-.
C*          (3,r,2) : (R*4) : m_(ij) + m_k  for r <--> (ij)
C*       Z(1-14)    : (R*8) : integration varialble in (0,1).
C* (Outputs)
C*       IPV(1,1) : (I*4) : a+
C*          (2,1) : (I*4) : b+
C*          (3,1) : (I*4) : c+
C*          (1,2) : (I*4) : a-
C*          (2,2) : (I*4) : b-
C*          (3,2) : (I*4) : c-
C*       XPH(1,1) : (R*4) : s.
C*          (2,1) : (R*4) : cos(theta_X-).
C*          (3,1) : (R*4) : phi_X-.
C*       XPH(1,2) : (R*4) : invariant mass squared for X+.
C*          (2,2) : (R*4) : cos(theta_a+) in X+ rest frame.
C*          (3,2) : (R*4) : phi_a+ in X+ rest frame.
C*       XPH(1,3) : (R*4) : invariant mass squared for a+b+.
C*          (2,3) : (R*4) : invariant mass squared for a+c+.
C*          (3,3) : (R*4) : phi_b+ in X+ rest frame.
C*       XPH(1,4) : (R*4) : invariant mass squared for X-.
C*          (2,4) : (R*4) : cos(theta_a-) in X- rest frame.
C*          (3,4) : (R*4) : phi_a- in X- rest frame.
C*       XPH(1,5) : (R*4) : invariant mass squared for a-b-.
C*          (2,5) : (R*4) : invariant mass squared for a-c-.
C*          (3,5) : (R*4) : phi_b- in X- rest frame..
C*       WAT      : (R*8) : Jacobian weight.
C* (Relation)
C*    Invokes UHQ2BW and UHQIJ3 in ttz_lib.
C* (Update Record)
C*   95/04/07  K.Fujii	Original version.
C*   95/04/09  K.Fujii	Restrict the Q3 and Q6 ranges to avioid
C*			extremely off-shell charginos, since such
C*			a case requires inclusions of diagrams 
C*			other than chargino pair productions.
C*   95/04/12  K.Fujii	Switch to the zero-width approximation
C*			when the width is below GMXMIN.
C*
CC**********************************************************************

      SUBROUTINE GETXPH(RS,AM,AMX,GMX,AMR,Z,IPV,XPH,WAT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      INTEGER*4  IPV(3,2)
      REAL   *4  RS, AM(*), AMX(2), GMX(2), AMR(3,3,2), XPH(3,*)
      REAL   *8  Z(*), WAT
C--
C  Local variables.
C--
      INTEGER*4  JPV(3,3,3)
      REAL   *8  WT
      DATA EPS    / 5.D0  /
      DATA GMXMIN / 1.D-2 /
      DATA NCALL  /    0  /
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
C--
         xPI   = ACOS(-1.D0)
         x2PI  = 2*xPI
         ENGM  = 10
C--
         JPV(1,1,2) = 3
         JPV(2,1,2) = 4
         JPV(3,1,2) = 5
         JPV(1,1,3) = 4
         JPV(2,1,3) = 3
         JPV(3,1,3) = 5
         JPV(1,2,1) = 3
         JPV(2,2,1) = 5
         JPV(3,2,1) = 4
         JPV(1,2,3) = 5
         JPV(2,2,3) = 3
         JPV(3,2,3) = 4
         JPV(1,3,1) = 4
         JPV(2,3,1) = 5
         JPV(3,3,1) = 3
         JPV(1,3,2) = 5
         JPV(2,3,2) = 4
         JPV(3,3,2) = 3
      ENDIF
C--
C  Set independent variables.
C     Z( 1) : m(X+)**2
C      ( 2) : m(X-)**2
C      ( 3) : m(+ab)**2
C      ( 4) : m(+ac)**2
C      ( 5) : m(-ab)**2
C      ( 6) : m(-ac)**2
C      ( 7) : cos(theta_X-)
C      ( 8) : phi_X-
C      ( 9) : cos(theta_a)     in X+ rest frame
C      (10) : phi_a            in X+ rest frame
C      (11) : phi_b            in X+ rest frame
C      (12) : cos(theta_a)     in X- rest frame
C      (13) : phi_a            in X- rest frame
C      (14) : phi_b            in X- rest frame
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Set invariant masses squared for X+.
C--
      AMSX  = AMX(1)
      GAMX  = GMX(1)
C--
      IF ( GAMX.GT.GMXMIN ) THEN
         Q3MN  = AM(3) + AM(4) + AM(5)
         Q3MX  = RS - (AM(6)+AM(7)+AM(8))
         IF ( AMSX.GT.Q3MN .AND. AMSX.LT.Q3MX ) THEN
            Q3MN  = MAX(Q3MN,AMSX-ENGM*GAMX)
            Q3MX  = MIN(Q3MX,AMSX+ENGM*GAMX)
         ENDIF
         CALL UHQ2BW(AMSX,GAMX,Q3MN,Q3MX,Z(5),Q32,WT)
         Q3    = SQRT(Q32)
      ELSE
         Q32 = AMSX*AMSX
         Q3  = AMSX
         WT  = xPI*AMSX*GAMX
      ENDIF      
      WAT   = WAT*WT
C--
C  And for X-.
C--
      AMSX  = AMX(2)
      GAMX  = GMX(2)
C--
      IF ( GAMX.GT.GMXMIN ) THEN
         Q6MN  = AM(6) + AM(7) + AM(8)
         Q6MX  = RS - Q3
         IF ( AMSX.GT.Q6MN .AND. AMSX.LT.Q6MX ) THEN
            Q6MN  = MAX(Q6MN,AMSX-ENGM*GAMX)
            Q6MX  = MIN(Q6MX,AMSX+ENGM*GAMX)
         ENDIF
         CALL UHQ2BW(AMSX,GAMX,Q6MN,Q6MX,Z(6),Q62,WT)
         Q6    = SQRT(Q62)
      ELSE
         Q62 = AMSX*AMSX
         Q6  = AMSX
         WT  = xPI*AMSX*GAMX
      ENDIF      
      WAT   = WAT*WT
C--
C  Then set invariant masses squared 
C  for daughter particles from X+.
C--
      CALL UHQIJ3(AM(1),JPV,AMR(1,1,1),EPS,Q3,Z(1),IPV(1,1),Q42,Q52,WT)
      WAT = WAT*WT
C--
C  And for daughter particles from X-.
C--
      CALL UHQIJ3(AM(4),JPV,AMR(1,1,2),EPS,Q6,Z(3),IPV(1,2),Q72,Q82,WT)
      WAT = WAT*WT
C--
      IPV(1,2) = IPV(1,2) + 3
      IPV(2,2) = IPV(2,2) + 3
      IPV(3,2) = IPV(3,2) + 3
C>>>
CCDBG       print *, ' IPV(*,1) = ', (IPV(K,1),K=1,3)
CCDBG       print *, '    (*,2) = ', (IPV(K,2),K=1,3)      
C>>>       
C--
C  Add angular variables.
C--
      XPH(1,1) = RS*RS
      XPH(2,1) = -1 + 2*Z(7)
      XPH(3,1) = x2PI*Z(8)
      XPH(1,2) = Q32
      XPH(2,2) = -1 + 2*Z(9)
      XPH(3,2) = x2PI*Z(10)
      XPH(1,3) = Q42
      XPH(2,3) = Q52
      XPH(3,3) = x2PI*Z(11)
      XPH(1,4) = Q62
      XPH(2,4) = -1 + 2*Z(12)
      XPH(3,4) = x2PI*Z(13)
      XPH(1,5) = Q72
      XPH(2,5) = Q82
      XPH(3,5) = x2PI*Z(14)
C>>>
CCDBG       PRINT *, ' Q32 = ', Q32, ' Q42 = ', Q42, ' Q52 = ', Q52
CCDBG       PRINT *, ' Q62 = ', Q62, ' Q72 = ', Q72, ' Q82 = ', Q82
C>>>
C--
C  Calculate Jacobian.
C--
      WAT = WAT*2**3*(x2PI)**5
C--
C  That's it.
C--
      RETURN
      END
