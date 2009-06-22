C*  09/06/18  K.Fujii           SPEVNT for ZZZ.
C*
C*
 
      SUBROUTINE SPEVNT(IRET)
 
      IMPLICIT    REAL*4 ( A-H, O-Z )
C--
C  Arguments.
C--
      INTEGER*4   IRET
C--
C  Commons.
C--
#include "ZZZCNF.inc"
#include "SMCONS.inc"
#include "SMPTAB.inc"
C--
C  Junk commons.
C--
      COMMON /SSWORK/ IBUF(10000)
      REAL   *4       RBUF(10000)
      EQUIVALENCE    (IBUF,RBUF)
C--
C  Local arrays.
C--
      REAL   *4       QV(0:3,10), QF(3,2,2)
C--
C  Initialize flags.
C--
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
C                 (14) = particle serial # of the parent daughter
C                 (15) =    not used
C                 (16) =    not used
C                 (17) = helicity
C                 (18) = colour single group ID
C                 (19) = parton-shower level flag
C                 (20) =    not used
C--
C     PVOUT(*,i) = 4-momentum of i-th parton
C
C                               /  3 f(1)
C                              /     _
C                             /___ 4 f(1)
C                            /
C      1 e-  \            Z1/         _
C             \  ######### /       /5 f(2)
C             ###############     /
C            #################___/
C            ################# Z2\
C             ###############     \
C             / ########## \       \6 f(2)
C      2 e+  /            Z \        _
C                            \___  7 f
C                             \
C                              \
C                               \  8 f
C--
      CALL PSUMxx(PVOUT(0,7),PVOUT(0,8),QV(0,2))
      CALL PSUMxx(PVOUT(0,3),PVOUT(0,4),QV(0,4))
      CALL PSUMxx(PVOUT(0,5),PVOUT(0,6),QV(0,5))
C--
      NPT = 0
C--
C  Store Z1.
C--
      CALL RDOTxx(QV(0,4),QV(0,4),VMW2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   23
      RBUF( 3) =   SQRT(VMW2)
      RBUF( 4) =   -1
      RBUF( 5) =   QV(1,4)
      RBUF( 6) =   QV(2,4)
      RBUF( 7) =   QV(3,4)
      RBUF( 8) =   QV(0,4)
      RBUF(12) =   2
      RBUF(13) =   4
      RBUF(14) =   0
      RBUF(18) =   1
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  Store Z2.
C--
      CALL RDOTxx(QV(0,5),QV(0,5),VMW2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =    23
      RBUF( 3) =   SQRT(VMW2)
      RBUF( 4) =    1
      RBUF( 5) =   QV(1,5)
      RBUF( 6) =   QV(2,5)
      RBUF( 7) =   QV(3,5)
      RBUF( 8) =   QV(0,5)
      RBUF(12) =   2
      RBUF(13) =   6
      RBUF(14) =   0
      RBUF(18) =   2
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  Store Z.
C--
      CALL RDOTxx(QV(0,2),QV(0,2),VMZ2)
      CALL UVZERO( 20, RBUF )
      NPT      =   NPT + 1
      RBUF( 1) =   NPT
      RBUF( 2) =   23
      RBUF( 3) =   SQRT(VMZ2)
      RBUF( 4) =   0
      RBUF( 5) =   QV(1,2)
      RBUF( 6) =   QV(2,2)
      RBUF( 7) =   QV(3,2)
      RBUF( 8) =   QV(0,2)
      RBUF(12) =   2
      RBUF(13) =   8
      RBUF(14) =   0
      RBUF(18) =   3
      CALL TBPUT(1,'Spring:Parton_List',NPT, 20, RBUF, IRT )
      IF ( NEVENT.LE.NPRINT ) THEN
         CALL SLPRNT(NEVENT,NPT,RBUF(1))
      ENDIF
C--
C  fermions from Z1.
C--
      NF = -1
      DO 200 IPT = 3, 4
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
         RBUF(14) =   1
         RBUF(17) =   IHLCMB(IPT)
         RBUF(18) =   1
         RBUF(19) =   (LQ-1)*101
         CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
         IF ( NEVENT.LE.NPRINT ) THEN
            CALL SLPRNT(NEVENT,NPT,RBUF(1))
         ENDIF
200   CONTINUE
C--
C  fermions from Z2.
C--
      NF = +1
      DO 300 IPT = 5, 6
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
         RBUF(14) =   2
         RBUF(17) =   IHLCMB(IPT)
         RBUF(18) =   2
         RBUF(19) =   (LQ-1)*201
         CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
         IF ( NEVENT.LE.NPRINT ) THEN
            CALL SLPRNT(NEVENT,NPT,RBUF(1))
         ENDIF
300   CONTINUE
C--
C  Fermions from Z.
C--
      NF = +1
      DO 400 IPT = 7, 8
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
         RBUF(14) =   3
         RBUF(17) =   IHLCMB(IPT)
         RBUF(18) =   3
         RBUF(19) =   (LQ-1)*301
         CALL TBPUT(1,'Spring:Parton_List',NPT,20,RBUF,IRT)
         IF ( NEVENT.LE.NPRINT ) THEN
            CALL SLPRNT(NEVENT,NPT,RBUF(1))
         ENDIF
400   CONTINUE
C--
C  That's it.
C--
      RETURN
      END
