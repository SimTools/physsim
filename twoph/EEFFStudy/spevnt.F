C* (Update Record)
C*  96/03/07  K.Fujii           FULL version for EEFF.
C*
C*
      SUBROUTINE SPEVNT

      IMPLICIT REAL*4 (A-H,O-Z)
C--
#include "USRPRM.inc"
#include "KINEM.inc"
#include "KINEM2.inc"
#include "KINEMS.inc"
#include "SMCONS.inc"
#include "SMPTAB.inc"
C--
      COMMON /SSWORK/ IBUF(10000)
      REAL   *4       RBUF(10000)
      EQUIVALENCE    (IBUF,RBUF)
C--
      INTEGER*4       IDPTCL(3,10), IHLCMB(10)
C--
      REAL   *4       QV(0:3,10), PV(0:3,10), QF(3,2,2)
      DATA NCALL      / 0 /
      DATA NPRINT     /10 /
      DATA NEVENT     / 0 /
C--
      REAL*8	      DRN      
C
C========< Entry Point >================================================
C
C--
C  (1) Reset IPET.
C--
      IPET = 0
C--
C  (2) Initialize constants.
C          IDPTCL(1,IPT) = IG  : generation
C                (2,IPT) = IT  : (1,2) = (up,down)
C                (3,IPT) = LQ  : (1,2) = (l,q)
C          QF(IG,IT,LQ)  = fermion charge
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         DO 10 IG = 1, 3
            QF(IG,1,1) =  0
            QF(IG,2,1) = -1
            QF(IG,1,2) =  2.D0/3
            QF(IG,2,2) = -1.D0/3
10       CONTINUE
      ENDIF
C--
      IDPTCL(1,1) = 1
      IDPTCL(2,1) = 2
      IDPTCL(3,1) = 1
      IDPTCL(1,2) = 1
      IDPTCL(2,2) = 2
      IDPTCL(3,2) = 1
      IDPTCL(1,3) = IGFR
      IDPTCL(2,3) = ITFR
      IDPTCL(3,3) = LQFR
      IDPTCL(1,4) = IGFR
      IDPTCL(2,4) = ITFR
      IDPTCL(3,4) = LQFR
C--
C  (3) Calculate kinematics.
C-- 
      COSTE1 = COSTH1/AK1
      COSTE2 = COSTH2/AK2
      SINTE1 = SINTH1/AK1
      SINTE2 = SINTH2/AK2
C--
C e-
C--
      P2      = E*AK2
      E2      = E*OMEGA2
C--
      QV(0,1) = E2
      QV(1,1) = P2*SINTE2*COSPH2
      QV(2,1) = P2*SINTE2*SINPH2
      QV(3,1) = P2*COSTE2
C--
C  e+
C--
      P1      = E*AK1
      E1      = E*OMEGA1
C--
      QV(0,2) = E1
      QV(1,2) = P1*SINTE1
      QV(2,2) = 0.0
      QV(3,2) = P1*COSTE1
C--
C  f
C--
      P3      = E*AKZ(IBX,IBCOS)
      COSTH   = COSTHZ(IBCOS)
      SINTH   = SINTHZ(IBCOS)
      E3      = E*OMEGAZ(IBX,IBCOS)
C--
      QV(0,3) = E3
      QV(1,3) = P3*SINTH*COSPHI
      QV(2,3) = P3*SINTH*SINPHI
      QV(3,3) = P3*COSTH
C--
C  fbar
C-- 
      QV(0,4) = - QV(0,1) - QV(0,2) - QV(0,3) + 2*E
      QV(1,4) = - QV(1,1) - QV(1,2) - QV(1,3)
      QV(2,4) = - QV(2,1) - QV(2,2) - QV(2,3)
      QV(3,4) = - QV(3,1) - QV(3,2) - QV(3,3)
C--
      P4      = SQRT( QV(1,4)**2 + QV(2,4)**2 + QV(3,4)**2 )
      E4      = QV(0,4)
C--
C  Generate an azimuthal angle to rotate the event.
C-- 
      PHIEVT  = 2*xPI*DRN(DUMM)
      CSPHI   = COS(PHIEVT)
      SNPHI   = SIN(PHIEVT)
      DO 20 IP = 1, 4
         PV(0,IP) = QV(0,IP)
         PV(1,IP) = QV(1,IP)*CSPHI - QV(2,IP)*SNPHI
         PV(2,IP) = QV(1,IP)*SNPHI + QV(2,IP)*CSPHI
         PV(3,IP) = QV(3,IP)
20    CONTINUE
C--
C  Generate fermion helicities.
C--
      IF ( DRN(DUMM).LE.0.5D0 ) THEN
         IHELF = -1
      ELSE
         IHELF = +1
      ENDIF
C--
      IHLCMB(1) =  0
      IHLCMB(2) =  0
      IHLCMB(3) =  IHELF
      IHLCMB(4) = -IHELF
C-- 
C  (4) Save events in TBS.
C--
C  TBPUT 'Spring:Header'.
C       Elm#1  IBUF(1) = Event #.
C                  (2) = Date
C                  (3) = Time
C                  (4) = not used
C                  (5) = E (e-)    MeV unit.
C                  (6) = E (e+)    MeV unit.
C                  (7) = Ebeam     MeV unit.
C                  (8) = Pol (e-)
C                  (9) = Pol (e+)
C                 (10) = not used.
C--
      NEVENT = NEVENT + 1
      IBUF(1) = NEVENT
      CALL UIDATE( IBUF(2) )
      IBUF(3) = TIME()
      IBUF(5) = E*1.E3
      IBUF(6) = E*1.E3
      IBUF(7) = E*1.E3
      IBUF(8) = POLE
      IBUF(9) = 0
      CALL TBPUT(1, 'Spring:Header', 1, 10, IBUF, IRT )
C--
C  TBPUT 'Spring:Parton_List'.
C       Elm#i  RBUF(1) = Particle serial number.
C                  (2) = Particle ID
C                  (3) = Mass (GeV)
C                  (4) = Charge
C                  (5) = Px(GeV)
C                  (6) = Py(GeV)
C                  (7) = Pz(GeV)
C                  (8) = E(GeV)
C                  (9) =    not used
C                 (10) =    not used
C                 (11) =    not used
C                 (12) = # daughter parton
C                 (13) = particle serial # of the 1st daughter
C                 (14) = particle serial # of the parent
C                 (15) =    not used
C                 (16) =    not used
C                 (17) = helicity
C                 (18) = colour single group ID
C                 (19) =    not used
C                 (20) =    not used
C--
C     PV(*,i) = 4-momentum of i-th parton
C
C        e-    -------+----------  1 e-
C                      \
C                       \
C                        +------   3  f
C                        |
C                        |            _
C                        +------   4  f
C                       /
C                      /
C        e+    -------+----------  2 e+
C
C--
      NPT = 0
C--
C  (a) Store final-state e+ and e-.
C--
      NF = -1
      DO 300 IPT = 1, 4
         CALL IDK2PD(IDPTCL(1,IPT),IDPDG)
         IG = IDPTCL(1,IPT)
         IT = IDPTCL(2,IPT)
         LQ = IDPTCL(3,IPT)
         NF = -NF
         CALL UVZERO( 20, RBUF )
         NPT      =   NPT + 1
         RBUF( 1) =   NPT
         RBUF( 2) =   NF*IDPDG
         RBUF( 3) =   AMF(IG,IT,LQ)
         RBUF( 4) =   NF*QF(IG,IT,LQ)
         RBUF( 5) =   PV(1,IPT)
         RBUF( 6) =   PV(2,IPT)
         RBUF( 7) =   PV(3,IPT)
         RBUF( 8) =   PV(0,IPT)
         RBUF(12) =   0
         RBUF(13) =   0
         RBUF(14) =   0
         RBUF(17) =   IHLCMB(IPT)
         RBUF(18) =   1
         RBUF(19) =   (LQ-1)*101
         CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
         IF ( NEVENT.LE.NPRINT ) THEN
            CALL SLPRNT(NEVENT,NPT,RBUF(1))
         ENDIF
300   CONTINUE
C--
C  That's it.
C--
      RETURN
      END

