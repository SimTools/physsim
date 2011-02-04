C*  2010/10/07  K.Fujii           SPEVNT for XNXN.
C*
 
      SUBROUTINE SPEVNT(IRET)
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
      PARAMETER        ( IxPRC = 1 )
#include "XNXNCN.inc"
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
      COMMON /SSWORK/ IBUF(10000)
      REAL   *4       RBUF(10000)
      EQUIVALENCE    (IBUF,RBUF)
      REAL   *8       QV(0:3,10), QF(3,2,2)
      DATA NCALL      / 0 /
      DATA NPRINT     /10 /
      DATA NEVENT     / 0 /
C
C========< Entry Point >================================================
C
C--
C  Reset IRET.
C--
      IRET = 0
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         DO 10 IG = 1, 3
            QF(IG,1,1) =  0.D0
            QF(IG,2,1) = -1.D0
            QF(IG,1,2) =  2.D0/3.D0
            QF(IG,2,2) = -1.D0/3.D0
10       CONTINUE
      ENDIF
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
      IBUF(5) = BMENGY(3)*1.D3
      IBUF(6) = BMENGY(4)*1.D3
      IBUF(7) = BMENGY(0)*1.D3
      IBUF(8) = POLE
      IBUF(9) = 0
      CALL TBPUT(1, 'Spring:Header', 1, 10, IBUF, IRT )
C--
C  TBPUT 'Spring:Header'.
C       Elm#2  IBUF( 1-18) = IDPTCL(3,6)
C                  (19-24) = IHLCMB
C              RBUF(25-33) = XPHASE(3,3)
C                  (   34) = DSGMDX
C                  (   35) = DPDEBM
C                  (   36) = nominal beam energy
C                  (   37) = Gaussian smeared E(e-)
C                  (   38) = Gaussian smeared E(e+)
C                  (   39) = E(e-) after beamstrahlung
C                  (   40) = E(e+) after beamstrahlung
C--
      CALL UVCOPY(18,IDPTCL,IBUF(1))
      CALL UVCOPY(6,IHLCMB(1),IBUF(19))
      CALL UVCOPY(9,XPHASE,RBUF(25))
      RBUF(34) = DSGMDX
      RBUF(35) = DPDEBM
      CALL UVCOPY(5,BMENGY(0),RBUF(36))
      CALL TBPUT(1, 'Spring:Header', 2, 40, IBUF, IRT )
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
C     PVOUT(*,i) = 4-momentum of i-th parton
C
C                               /  3 tau-/+
C                              /
C                             /  
C                            /\  
C      1 e-  \          X0i /  \
C             \            /    \  
C              \          /      \ 4 stau+/-
C               \________/
C               /        \
C              /          \      / 6 stau-/+
C             /            \    /
C      2 e+  /          X0j \  /
C                            \/
C                             \
C                              \
C                               \  5 tau+/-
C--
      NPT = 0
C--
C  (1) Store X0_i.
C--
      CALL PSUMxx(PVOUT(0,3),PVOUT(0,4),QV(0,3))
      CALL RDOTxx(QV(0,3),QV(0,3),VMT2)
C--
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   222000
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =   0.D0
      RBUF( 5) =   QV(1,3)
      RBUF( 6) =   QV(2,3)
      RBUF( 7) =   QV(3,3)
      RBUF( 8) =   QV(0,3)
      RBUF(12) =   2
      RBUF(13) =   3
      RBUF(14) =   0
      RBUF(18) =   1
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (2) Store X0_j.
C--
      CALL PSUMxx(PVOUT(0,5),PVOUT(0,6),QV(0,4))
      CALL RDOTxx(QV(0,4),QV(0,4),VMT2)
C--
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   222000
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =   0.D0
      RBUF( 5) =   QV(1,4)
      RBUF( 6) =   QV(2,4)
      RBUF( 7) =   QV(3,4)
      RBUF( 8) =   QV(0,4)
      RBUF(12) =   2
      RBUF(13) =   5
      RBUF(14) =   0
      RBUF(18) =   1
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (3) matter fermion from X0_i.
C--
      IPT = 3
      CALL IDK2PD(IDPTCL(1,IPT),IDPDG)
      IG = IDPTCL(1,IPT)
      IT = IDPTCL(2,IPT)
      LQ = IDPTCL(3,IPT)
      IF (IDPTCL(2,4).LT.0) THEN
         NF = -1
      ELSE
         NF = +1
      ENDIF
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   NF*IDPDG
      RBUF( 3) =   AMF(IG,IT,LQ)
      RBUF( 4) =   NF*QF(IG,IT,LQ)
      RBUF( 5) =   PVOUT(1,IPT)
      RBUF( 6) =   PVOUT(2,IPT)
      RBUF( 7) =   PVOUT(3,IPT)
      RBUF( 8) =   PVOUT(0,IPT)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   1
      RBUF(17) =   IHLCMB(IPT)
      RBUF(18) =   3
      RBUF(19) =   (LQ-1)*101
      CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (4) sfermion from X0_i.
C--
      IPT = 4
#if 1
      IDPDG = 1000015
#else
      IDPDG = 13
#endif
      IG = IDPTCL(1,IPT)
      IT = IDPTCL(2,IPT)
      LQ = IDPTCL(3,IPT)
      IF (IDPTCL(2,4).LT.0) THEN
         NF = -1
      ELSE
         NF = +1
      ENDIF
      ITA = IABS(IT)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   NF*IDPDG
      RBUF( 3) =   SFM3M(2)
      RBUF( 4) =   NF*QF(IG,ITA,LQ)
      RBUF( 5) =   PVOUT(1,IPT)
      RBUF( 6) =   PVOUT(2,IPT)
      RBUF( 7) =   PVOUT(3,IPT)
      RBUF( 8) =   PVOUT(0,IPT)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   1
      RBUF(17) =   IHLCMB(IPT)
      RBUF(18) =   3
      RBUF(19) =   (LQ-1)*101
      CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (5) matter fermion from X0_i.
C--
      IPT = 5
      CALL IDK2PD(IDPTCL(1,IPT),IDPDG)
      IG = IDPTCL(1,IPT)
      IT = IDPTCL(2,IPT)
      LQ = IDPTCL(3,IPT)
      IF (IDPTCL(2,6).LT.0) THEN
         NF = +1
      ELSE
         NF = -1
      ENDIF
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   NF*IDPDG
      RBUF( 3) =   AMF(IG,IT,LQ)
      RBUF( 4) =   NF*QF(IG,IT,LQ)
      RBUF( 5) =   PVOUT(1,IPT)
      RBUF( 6) =   PVOUT(2,IPT)
      RBUF( 7) =   PVOUT(3,IPT)
      RBUF( 8) =   PVOUT(0,IPT)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   2
      RBUF(17) =   IHLCMB(IPT)
      RBUF(18) =   4
      RBUF(19) =   (LQ-1)*101
      CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (6) sfermion from X0_i.
C--
      IPT = 6
#if 1
      IDPDG = 1000015
#else
      IDPDG = 13
#endif
      IG = IDPTCL(1,IPT)
      IT = IDPTCL(2,IPT)
      LQ = IDPTCL(3,IPT)
      IF (IDPTCL(2,6).LT.0) THEN
         NF = +1
      ELSE
         NF = -1
      ENDIF
      ITA = IABS(IT)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   NF*IDPDG
      RBUF( 3) =   SFM3M(2)
      RBUF( 4) =   NF*QF(IG,ITA,LQ)
      RBUF( 5) =   PVOUT(1,IPT)
      RBUF( 6) =   PVOUT(2,IPT)
      RBUF( 7) =   PVOUT(3,IPT)
      RBUF( 8) =   PVOUT(0,IPT)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   2
      RBUF(17) =   IHLCMB(IPT)
      RBUF(18) =   4
      RBUF(19) =   (LQ-1)*101
      CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
