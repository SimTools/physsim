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
C*       RS         : (R*8) : sqrt(s).
C*       AM(i)      : (R*8) : mass  of i.
C*       AMX(1)     : (R*8) : mass  of X0i.
C*          (2)     : (R*8) : mass  of X0j.
C*       GMX(1)     : (R*8) : width of X0i.
C*          (2)     : (R*8) : width of X0j.
C*       AMR(1,r)   : (R*8) : mass of resonance r        
C*          (2,r)   : (R*8) : width of resonance r           for X0j.
C*          (3,r)   : (R*8) : m_(ij) + m_k  for r <--> (ij)
C*       Z(1-8)     : (R*8) : integration varialble in (0,1).
C* (Outputs)
C*       IPV(1)   : (I*4) : a
C*          (2)   : (I*4) : b
C*          (3)   : (I*4) : c
C*       XPH(1,1) : (R*8) : s.
C*          (2,1) : (R*8) : cos(theta_X0j).
C*          (3,1) : (R*8) : phi_X0j.
C*       XPH(1,2) : (R*8) : invariant mass squared for X0j.
C*          (2,2) : (R*8) : cos(theta_a) in X0j rest frame.
C*          (3,2) : (R*8) : phi_a in X0j rest frame.
C*       XPH(1,3) : (R*8) : invariant mass squared for ab.
C*          (2,3) : (R*8) : invariant mass squared for ac.
C*          (3,3) : (R*8) : phi_b in X0j rest frame.
C*       WAT      : (R*8) : Jacobian weight.
C* (Relation)
C*    Invokes UHQ2BW and UHQIJ3.
C* (Update Record)
C*   03/08/28  K.Fujii	Original version.
C*
CC**********************************************************************

      SUBROUTINE GETXPH(RS,AM,AMX,GMX,AMR,Z,IPV,XPH,WAT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      INTEGER*4  IPV(3)
      REAL   *8  RS, AM(*), AMX(2), GMX(2), AMR(3,3), XPH(3,*)
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
C         ENGM  = 1000
         ENGM  = 100
C--
C  Resonace index and table
C        JPV(1,r_1,r_2) : (I*4) : i    
C           (2,r_1,r_2) : (I*4) : j   r_1 <--> (ij)
C           (3,r_1,r_2) : (I*4) : k   r_2 <--> (ik)
C        AMR(1,r)       : (R*8) : mass of resonance r
C           (2,r)       : (R*8) : width of resonance r
C           (3,r)       : (R*8) : m_(ij) + m_k  for r <--> (ij)
C--
         JPV(1,1,2) = 6
         JPV(2,1,2) = 4
         JPV(3,1,2) = 5
         JPV(1,1,3) = 4
         JPV(2,1,3) = 6
         JPV(3,1,3) = 5
         JPV(1,2,1) = 6
         JPV(2,2,1) = 5
         JPV(3,2,1) = 4
         JPV(1,2,3) = 5
         JPV(2,2,3) = 6
         JPV(3,2,3) = 4
         JPV(1,3,1) = 4
         JPV(2,3,1) = 5
         JPV(3,3,1) = 6
         JPV(1,3,2) = 5
         JPV(2,3,2) = 4
         JPV(3,3,2) = 6
      ENDIF
C--
C  Set independent variables.
C     Z( 1) : m(ab)**2
C      ( 2) : m(ac)**2
C      ( 3) : m(X0j)**2
C      ( 4) : cos(theta_X0j)
C      ( 5) : phi_X0j
C      ( 6) : cos(theta_a)     in X0j rest frame
C      ( 7) : phi_a            in X0j rest frame
C      ( 8) : phi_b            in X0j rest frame
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Set invariant masses squared for X0j.
C--
      JX = 2
C--
      AMSX  = AMX(JX)
      GAMX  = GMX(JX)
C--
      IF ( GAMX.GT.GMXMIN ) THEN
         Q3MN  = AM(4) + AM(5) + AM(6)
         Q3MX  = RS - AM(3)
         IF ( AMSX.GT.Q3MN .AND. AMSX.LT.Q3MX ) THEN
            Q3MN  = MAX(Q3MN,AMSX-ENGM*GAMX)
            Q3MX  = MIN(Q3MX,AMSX+ENGM*GAMX)
         ENDIF
         CALL UHQ2BW(AMSX,GAMX,Q3MN,Q3MX,Z(3),Q32,WT)
         Q3    = SQRT(Q32)
      ELSE
         Q32 = AMSX*AMSX
         Q3  = AMSX
         WT  = xPI*AMSX*GAMX
      ENDIF      
      WAT   = WAT*WT
C--
C  Then set invariant masses squared 
C  for daughter particles from X0j.
C--
      CALL UHQIJ3(AM,JPV,AMR,EPS,Q3,Z(1),IPV,Q42,Q52,WT)
      WAT = WAT*WT
C>>>
CCCDBG       print *, ' IPV(*) = ', (IPV(K),K=1,3)
C>>>       
C--
C  Add angular variables.
C--
      XPH(1,1) = RS*RS
      XPH(2,1) = -1 + 2*Z(4)
      XPH(3,1) = x2PI*Z(5)
      XPH(1,2) = Q32
      XPH(2,2) = -1 + 2*Z(6)
      XPH(3,2) = x2PI*Z(7)
      XPH(1,3) = Q42
      XPH(2,3) = Q52
      XPH(3,3) = x2PI*Z(8)
C>>>
CCCDBG       PRINT *, ' Q32 = ', Q32, ' Q42 = ', Q42, ' Q52 = ', Q52
C>>>
C--
C  Calculate Jacobian.
C--
      WAT = WAT * 2**2 * (x2PI)**3
C--
C  That's it.
C--
      RETURN
      END
