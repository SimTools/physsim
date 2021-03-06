CC**********************************************************************
C*
C*===================----==
C* Subroutine SPEVNT(IRET)
C*===================----==
C* (Purpose)
C*   To enerate One event.
C* (Output)
C*   IRET :(I*4): = 0     Normal return.
C*                > 0     Warning, continue processing.
C*                < 0     Error, skip the rest.
C*                > -999  terminates the job imediately.
C* (Update Record)
C*  2010/11/29  K.Fujii          Derived from TTStudy/spevnt.F
C*
CC**********************************************************************
 
      SUBROUTINE SPEVNT(IRET)
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
      PARAMETER        ( IxPRC = 1 )
      INCLUDE 'TBWCNF.inc'
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "BEMCNS.inc"
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
C                                    _
C                               /  3 b
C                              /
C                             /      _
C                         _  /\W-/ 4 f(up)
C      1 e-  \            t /  \/
C             \            /    \
C              \          /      \ 5 f(down)
C               \________/
C               /        \           _
C              /          \      / 8 f(down)
C             /            \    /
C      2 e+  /             t\  /\
C                            \/W+\ 7 f(up)
C                             \
C                              \
C                               \  6 b(down)
C--
      CALL PSUMxx(PVOUT(0,4),PVOUT(0,5),QV(0,7))
      CALL PSUMxx(QV(0,7),PVOUT(0,3),QV(0,3))
      CALL PSUMxx(PVOUT(0,7),PVOUT(0,8),QV(0,8))
      CALL PSUMxx(QV(0,8),PVOUT(0,6),QV(0,4))
C--
      NPT = 0
C--          _
C  (1) Store t.
C--
      CALL RDOTxx(QV(0,3),QV(0,3),VMT2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   -6
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =   -2.D0/3.D0
      RBUF( 5) =   QV(1,3)
      RBUF( 6) =   QV(2,3)
      RBUF( 7) =   QV(3,3)
      RBUF( 8) =   QV(0,3)
      RBUF(12) =   2
      RBUF(13) =   3
      RBUF(14) =   0
      RBUF(18) =   1
      RBUF(19) =   0
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (2) Store t.
C--
      CALL RDOTxx(QV(0,4),QV(0,4),VMT2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =    6
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =    2.D0/3.D0
      RBUF( 5) =   QV(1,4)
      RBUF( 6) =   QV(2,4)
      RBUF( 7) =   QV(3,4)
      RBUF( 8) =   QV(0,4)
      RBUF(12) =   2
      RBUF(13) =   5
      RBUF(14) =   0
      RBUF(18) =   1
      RBUF(19) =   0
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--          _      _
C  (3) Store q from t.
C--
      CALL IDK2PD(IDPTCL(1,3),IDPDG)
      IG = IDPTCL(1,3)
      IT = IDPTCL(2,3)
      LQ = IDPTCL(3,3)
      NF = -1
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   NF*IDPDG
      RBUF( 3) =   AMF(IG,IT,LQ)
      RBUF( 4) =   NF*QF(IG,IT,LQ)
      RBUF( 5) =   PVOUT(1,3)
      RBUF( 6) =   PVOUT(2,3)
      RBUF( 7) =   PVOUT(3,3)
      RBUF( 8) =   PVOUT(0,3)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   1
      RBUF(17) =   IHLCMB(3)
      RBUF(18) =   2
      RBUF(19) =   101
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--                  _
C  (4) Store W- from t.
C--
      CALL RDOTxx(QV(0,7),QV(0,7),VMT2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =  -24
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =  -1.D0
      RBUF( 5) =   QV(1,7)
      RBUF( 6) =   QV(2,7)
      RBUF( 7) =   QV(3,7)
      RBUF( 8) =   QV(0,7)
      RBUF(12) =   2
      RBUF(13) =   7
      RBUF(14) =   1
      RBUF(18) =   3
      RBUF(19) =   101
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (5) Store q from t.
C--
      CALL IDK2PD(IDPTCL(1,6),IDPDG)
      IG = IDPTCL(1,6)
      IT = IDPTCL(2,6)
      LQ = IDPTCL(3,6)
      NF = +1
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   NF*IDPDG
      RBUF( 3) =   AMF(IG,IT,LQ)
      RBUF( 4) =   NF*QF(IG,IT,LQ)
      RBUF( 5) =   PVOUT(1,6)
      RBUF( 6) =   PVOUT(2,6)
      RBUF( 7) =   PVOUT(3,6)
      RBUF( 8) =   PVOUT(0,6)
      RBUF(12) =   0
      RBUF(13) =   0
      RBUF(14) =   2
      RBUF(17) =   IHLCMB(6)
      RBUF(18) =   2
      RBUF(19) =   102
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (6) Store W+ from t.
C--
      CALL RDOTxx(QV(0,8),QV(0,8),VMT2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   24
      RBUF( 3) =   SQRT(VMT2)
      RBUF( 4) =   1.D0
      RBUF( 5) =   QV(1,8)
      RBUF( 6) =   QV(2,8)
      RBUF( 7) =   QV(3,8)
      RBUF( 8) =   QV(0,8)
      RBUF(12) =   2
      RBUF(13) =   9
      RBUF(14) =   2
      RBUF(18) =   4
      RBUF(19) =   102
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  (7-8) fermions from W-.
C--
      NF = 1
      DO 200 IPT = 4, 5
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
         RBUF( 5) =   PVOUT(1,IPT)
         RBUF( 6) =   PVOUT(2,IPT)
         RBUF( 7) =   PVOUT(3,IPT)
         RBUF( 8) =   PVOUT(0,IPT)
         RBUF(12) =   0
         RBUF(13) =   0
         RBUF(14) =   4
         RBUF(17) =   IHLCMB(IPT)
         RBUF(18) =   3
         RBUF(19) =   (LQ-1)*201
         CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
         IF ( NEVENT.LE.NPRINT ) THEN
            CALL SLPRNT(NEVENT,NPT,RBUF(1))
         ENDIF
200   CONTINUE
C--
C  (9-10) Fermions from W+.
C--
      NF = -1
      DO 300 IPT = 7, 8
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
         RBUF( 5) =   PVOUT(1,IPT)
         RBUF( 6) =   PVOUT(2,IPT)
         RBUF( 7) =   PVOUT(3,IPT)
         RBUF( 8) =   PVOUT(0,IPT)
         RBUF(12) =   0
         RBUF(13) =   0
         RBUF(14) =   6
         RBUF(17) =   IHLCMB(IPT)
         RBUF(18) =   4
         RBUF(19) =   (LQ-1)*301
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
