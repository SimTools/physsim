C* (Update Record)
C*  92/06/21  K.Fujii           SPEVNT for SFSF.
C*  95/02/12  K.Fujii		Now stores lepton helicities.
C*  95/03/21  K.Fujii           For BASES V5.1.
C*  95/05/05  K.Fujii           Special treatment for nonzero IDOTAU.
C*
 
      SUBROUTINE SPEVNT(IRET)
 
      IMPLICIT    REAL*4 ( A-H, O-Z )
      PARAMETER        ( IxPRC = 1 )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SFSFCN.inc"
      COMMON /SSWORK/ IBUF(10000)
      REAL   *4       RBUF(10000)
      EQUIVALENCE    (IBUF,RBUF)
      REAL   *4       QV(0:3,10), QF(3,2,2)
      DATA NCALL      / 0 /
      DATA NPRINT     /10 /
      DATA NEVENT     / 0 /
C--
C  Particle ID for sfermions and LSP.
C--
      INTEGER*4       IDSF(2,3), IDX01
      DATA IDSF       / 1000011, 2000011,
     .                  1000013, 2000013,
     .                  1000015, 2000015 / 
      DATA IDX01      / 1000022 /
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
            QF(IG,1,1) =  0.
            QF(IG,2,1) = -1.
            QF(IG,1,2) =  2./3.
            QF(IG,2,2) = -1./3.
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
      IBUF(5) = BMENGY(3)*1.E3
      IBUF(6) = BMENGY(4)*1.E3
      IBUF(7) = BMENGY(0)*1.E3
      IBUF(8) = POLE
      IBUF(9) = 0
      CALL TBPUT(1, 'Spring:Header', 1, 10, IBUF, IRT )
C--
C  TBPUT 'Spring:Header'.
C       Elm#2  IBUF( 1-24) = IDPTCL(3,8)
C                  (25-32) = IHLCMB
C              RBUF(33-47) = XPHASE(3,5)
C                  (   48) = DSGMDX
C                  (   49) = DPDEBM
C                  (   50) = nominal beam energy
C                  (   51) = Gaussian smeared E(e-)
C                  (   52) = Gaussian smeared E(e+)
C                  (   53) = E(e-) after beamstrahlung
C                  (   54) = E(e+) after beamstrahlung
C--
      CALL UVCOPY(24,IDPTCL,IBUF(1))
      CALL UVCOPY(8,IHLCMB(1),IBUF(25))
      CALL UVCOPY(15,XPHASE,RBUF(33))
      RBUF(48) = DSGMDX
      RBUF(49) = DPDEBM
      CALL UVCOPY(5,BMENGY(0),RBUF(50))
      CALL TBPUT(1, 'Spring:Header', 2, 54, IBUF, IRT )
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
C                               /  3 X0
C                              /
C                             /
C                            /\
C      1 e-  \           S- /  \
C             \            /    \
C              \          /      \ 4 l-
C               \________/
C               /        \
C              /          \      / 6 l+
C             /            \    /
C      2 e+  /           S+ \  /
C                            \/
C                             \
C                              \
C                               \  5 X0
C--
      CALL PSUMxx(PVOUT(0,3),PVOUT(0,4),QV(0,3))
      CALL PSUMxx(PVOUT(0,5),PVOUT(0,6),QV(0,4))
C--
      NPT = 0
C--
C  (1) Store S-.
C--
      CALL RDOTxx(QV(0,3),QV(0,3),VMT2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   IDSF(IHANDM,IGENSF)
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =  -1
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
C  (2) Store S+.
C--
      CALL RDOTxx(QV(0,4),QV(0,4),VMT2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =  -IDSF(IHANDP,IGENSF)
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =  +1
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
C  (3) Store X0 from S-.
C--
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   IDX01
      RBUF( 3) =   SZM(1)
      RBUF( 4) =   0
      RBUF( 5) =   PVOUT(1,3)
      RBUF( 6) =   PVOUT(2,3)
      RBUF( 7) =   PVOUT(3,3)
      RBUF( 8) =   PVOUT(0,3)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   1
      RBUF(18) =   1
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (4) matter fermion from S-.
C--
      IPT = 4
      CALL IDK2PD(IDPTCL(1,IPT),IDPDG)
      IG = IDPTCL(1,IPT)
      IT = IDPTCL(2,IPT)
      LQ = IDPTCL(3,IPT)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   IDPDG
      RBUF( 3) =   AMF(IG,IT,LQ)
      RBUF( 4) =   QF(IG,IT,LQ)
      RBUF( 5) =   PVOUT(1,IPT)
      RBUF( 6) =   PVOUT(2,IPT)
      RBUF( 7) =   PVOUT(3,IPT)
      RBUF( 8) =   PVOUT(0,IPT)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   1
      IF ( IDOTAU.EQ.0 ) THEN
         RBUF(17) = IHLCMB(4)
      ELSE
         RBUF(17) = HTAUM
      ENDIF
      RBUF(18) =   2
      CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (5) Store X0 from S+.
C--
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   IDX01
      RBUF( 3) =   SZM(1)
      RBUF( 4) =   0
      RBUF( 5) =   PVOUT(1,5)
      RBUF( 6) =   PVOUT(2,5)
      RBUF( 7) =   PVOUT(3,5)
      RBUF( 8) =   PVOUT(0,5)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   2
      RBUF(18) =   3
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (6) matter fermion from S+.
C--
      IPT = 6
      CALL IDK2PD(IDPTCL(1,IPT),IDPDG)
      IG = IDPTCL(1,IPT)
      IT = IDPTCL(2,IPT)
      LQ = IDPTCL(3,IPT)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =  -IDPDG
      RBUF( 3) =   AMF(IG,IT,LQ)
      RBUF( 4) =  -QF(IG,IT,LQ)
      RBUF( 5) =   PVOUT(1,IPT)
      RBUF( 6) =   PVOUT(2,IPT)
      RBUF( 7) =   PVOUT(3,IPT)
      RBUF( 8) =   PVOUT(0,IPT)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   2
      IF ( IDOTAU.EQ.0 ) THEN
         RBUF(17) = IHLCMB(6)
      ELSE
         RBUF(17) = HTAUP
      ENDIF
      RBUF(18) =   4
C--
      CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
