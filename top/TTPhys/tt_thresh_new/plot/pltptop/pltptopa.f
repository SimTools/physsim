C  C*
C  C* (Update Record)
C  C*   92/03/19  K.Fujii       Parametrization of Cho's new potential
C  C*                           parameters.
C  C*   92/03/19  K.Fujii       This version uses the exact values
C  C*                           obtained by Cho.
C  C*
C
C        SUBROUTINE GTPPAR(X,PPAR)
C
C        IMPLICIT REAL*8 ( A-H, O-Z )
C        REAL*8   X, PPAR(3)
C  C--
C        PARAMETER      ( NPT = 5 )
C        REAL*4           YDAT(NPT,3), DYDAT(NPT,3)
C        DATA YDAT  / 0.20823, 0.23347, 0.23530, 0.20681, 0.16754,
C       .             3.9620 , 3.8077 , 3.7397 , 3.7453 , 3.7485 ,
C       .             0.35914, 0.35398, 0.35719, 0.37349, 0.40100 /
C        DATA DYDAT / 0.01028, 0.01095, 0.00954, 0.00572, 0.00203,
C       .             0.27724, 0.28706, 0.35451, 0.40334, 0.29333,
C       .             0.02136, 0.02130, 0.02466, 0.02769, 0.02078 /
C  C
C  C==================< Entry Point >===================================
C  C
C  C--
C  C  Calculate PPAR.
C  C--
C        IF ( X.GE.0.099D0 .AND. X.LE.0.101D0 ) THEN
C           J = 1
C        ELSE IF ( X.GE.0.109D0 .AND. X.LE.0.111D0 ) THEN
C           J = 2
C        ELSE IF ( X.GE.0.119D0 .AND. X.LE.0.121D0 ) THEN
C           J = 3
C        ELSE IF ( X.GE.0.129D0 .AND. X.LE.0.131D0 ) THEN
C           J = 4
C        ELSE IF ( X.GE.0.139D0 .AND. X.LE.0.141D0 ) THEN
C           J = 5
C        ELSE
C           PRINT *, ' ERROR in GTPPAR: ALFS = ', X
C        ENDIF
C  C>>>
C  C     PPAR(1) = YDAT(J,1) - 5*DYDAT(J,1)
C        PPAR(1) = YDAT(J,1)
C  C     PPAR(1) = YDAT(J,1) + 5*DYDAT(J,1)
C  C>>>
C        PPAR(2) = YDAT(J,2)
C        PPAR(3) = YDAT(J,3)
C  C--
C  C  That's it.
C  C--
C        RETURN
C        END
C
C=======================================================================
C  Subroutines necessary to run this job.
C=======================================================================
C
C--
C  QCD corrections.
C--
CINCLUDE GTCEFF
CINCLUDE GMTWMB
CINCLUDE QCDPOT
CINCLUDE ADJ1
CINCLUDE ADJ2
CINCLUDE QCDS1
CINCLUDE QCDS2
CINCLUDE QCDS3
CINCLUDE BNDCND
CINCLUDE RKUTAS
CINCLUDE QCDP1
CINCLUDE QCDP2
CINCLUDE QCDP3
CINCLUDE RKUTAP
CINCLUDE DELFN
CINCLUDE RUNWID
CINCLUDE LNT0FI
CINCLUDE LNT1FI
CINCLUDE NTGRD1
CINCLUDE NTGRD2
CINCLUDE NTGRD3
CINCLUDE LNT0WA
CINCLUDE PHS0
CINCLUDE LNT2WA
CINCLUDE PHS2
CINCLUDE TBWCORR
CINCLUDE FX
CINCLUDE GTQCDP
CINCLUDE GAMMA0
C--
C  Cross section and asymmetry.
C--
CINCLUDE DSGDPF
CINCLUDE DSGDPA
CINCLUDE DSGDP0
CINCLUDE GRQDSG
CINCLUDE GRQINT
CINCLUDE GRQADJ
CINCLUDE GRQPTP
C--
C  ISR and beam effects.
C       RSHDIS0 : ISR only.
C       RSHDIS  : ISR + Beam effects.
C--
CINCLUDE UDSRCH
CINCLUDE RSHDIS
CINCLUDE RSHDIS0
C--
C  E_peak hunting.
C--
CINCLUDE EPKHNT
CINCLUDE SGTTHR2
CINCLUDE RKUTTA2
CINCLUDE EI
CINCLUDE F1
CINCLUDE FH
CINCLUDE GTGAMT
CINCLUDE SETPRMC
C
C=======================================================================
C  Main programs.
C=======================================================================
C
C*
C* Main program to draw the dsig/dp for various ALFS.
C* (Update Record)
C*   94/03/30  K.Fujii       Full o(alpha_s) version.
C*                           This version plots ptop for given DE's.
C*
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      PARAMETER     ( MXxNP = 200, MXxNE = 20, MXxNX = 10 )
      REAL     *8   PTDATA(10,0:MXxNP,0:MXxNE,0:MXxNX)
C--
C  E1S is an approximate 1S peak position.
C--
      REAL     *8   E1S
      PARAMETER     ( NE = 2 )
      REAL     *8   DEDATA(0:NE)
C      PARAMETER     ( NE = 0 )
C      REAL     *8   DEDATA(0:NE)
C--
      CHARACTER*8   JOIN(0:3)
C>>>
      DATA E1S      / -2.5D0 /
      DATA DEDATA   /   0.D0, 2.D0, 4.D0 /
C     DATA DEDATA   /   2.D0 /
C>>>
C--
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOU      / 21 /
C
C========< Entry Point >================================================
C
C--
C  Set parameters.
C--
      ALPSMN = 0.11D0
      ALPSMX = 0.13D0
      NA     = 2
      DA     = (ALPSMX-ALPSMN)/NA
C>>>
      NA     = 0
      ALPSMN = 0.12D0
C>>>
      VTB2   = 1.D0
      AMST   = 170.D0
      AMSH   = 1.D8
      BETH   = 0
C--
      ALP    = 1/128.0D0
      SN2W   = 0.23D0
      AMSZ   = 91.17D0
      GAMZ   = 2.5D0
      AMSW   = 80.0D0
      AMSB   = 5.0D0
C--
C  Initialize radiation function.
C--
      CALL ISRINT(0,2*AMST)
C--
C  Momentum range.
C--
      NQ     = 50
      QMN    =  0
      QMX    = 50
      DQ     = (QMX-QMN)/NQ
C--
C  Allocate output data set.
C--
C     CALL UALCPS(LOU,'TKSF.@@.@PTOP@.ALFS.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@PTOP@.ALFS.BM.TDR','RENEW','FB',30,IRT)
C--
C  Loop over alpha_s.
C--
      MODE = 0
      DO 10000 IA = 0, NA
         ALPS  = ALPSMN + IA*DA
         PRINT *, ' alpha_s  = ', ALPS
C--
C  First find 1S peak position.
C--
         EPK = E1S
         CALL EPKHNT( 0, ALP,SN2W,AMSZ,AMSW,AMSB,
     .                AMST,ALPS,VTB2,AMSH,BETH, EPK,SGPK, NTRY )
         PRINT *, ' NTRY, EPK, SGPK =', NTRY, EPK, SGPK
C--
C  Loop over energy.
C--
         DO 1000 IE = 0, NE
            E   = EPK + DEDATA(IE)
            RS     = 2*AMST +  E
C--
C  |p| distribution.
C--
            SG    = 0
            DO 100 IQ = 0, NQ
               IF ( IQ.EQ.0 .OR. IQ.EQ.NQ ) THEN
                  WT = 1
               ELSE IF ( MOD(IQ,2).EQ.0 ) THEN
                  WT = 2
               ELSE
                  WT = 4
               ENDIF
               Q    = QMN + IQ*DQ
C--
               CALL DSGDPF(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,
     .                     AMST, VTB2,  E,Q,DSG)
               MODE = 4
               SG   = SG + DSG*DQ*WT/3
C>>>
C               print *, ' Q, DSG = ', Q, DSG
C>>>
C--
               PTDATA(1,IQ,IE,IA) = Q
               PTDATA(2,IQ,IE,IA) = DSG
               PTDATA(3,IQ,IE,IA) = E
               PTDATA(4,IQ,IE,IA) = AMST
               PTDATA(5,IQ,IE,IA) = ALPS
               PTDATA(6,IQ,IE,IA) = VTB2
               PTDATA(7,IQ,IE,IA) = AMSH
               PTDATA(8,IQ,IE,IA) = BETH
100         CONTINUE
            MODE  = 3
            PRINT *, ' sqrt(s), sig_tt = ', RS, SG, ' pb'
1000     CONTINUE
         MODE  = 1
10000 CONTINUE
C--
C  Prepare TOPDRAW files.
C--
      WRITE(LOU,'(''( =========='')')
      WRITE(LOU,'(''(  dsig/dp'')')
      WRITE(LOU,'(''( =========='')')
      WRITE(LOU,'(''('')')
C--
      DO 20000 IE = 0, NE
         E      = PTDATA(3,0,IE,0)
         DE     = DEDATA(IE)
C--
         WRITE(LOU,'(''( =================='')')
         WRITE(LOU,'(''( DE = '',1F8.4)') DE
         WRITE(LOU,'(''( =================='')')
         WRITE(LOU,'(''('')')
         WRITE(LOU,'(''NEWFRAME'')')
         WRITE(LOU,'(''SET FONT DUPLEX'')')
         WRITE(LOU,'(''SET INTENSITY 4'')')
         WRITE(LOU,'(''('')')
         WRITE(LOU,'(''SET WINDOW X  2.5 12.5 Y 1.8  9.2'')')
         WRITE(LOU,'(''SET LIMIT X '',2E15.8)') QMN, QMX
         WRITE(LOU,'(''SET LIMIT Y 0. 0.03'')')
         WRITE(LOU,'(''(SET SCALE Y LOG'')')
         WRITE(LOU,'(''SET LABELS SIZE 3.0'')')
C--
         IF ( IE.EQ.0 ) THEN
            WRITE(LOU,'(''TITLE 3.0 8.7 SIZE 3 ''''*) E = E01S1'')')
            WRITE(LOU,'(''CASE                 ''''        X  X'')')
            WRITE(LOU,'(''TITLE 5.2 8.7 SIZE 3 ''''+'',1F5.2,''GeV'')')
     .                     DE
         ELSE
            WRITE(LOU,'(''TITLE 3.0 8.7 SIZE 3 ''''*) E ='',1F4.1,
     .                  '' GeV'')') E
         ENDIF
         WRITE(LOU,'(''TITLE 0.5 3.2 ANGLE 90 SIZE 4 ''''dS0tt1/d|p|'',
     .            '' (pb/GeV)'')')
         WRITE(LOU,'(''CASE                       '''' GX  X     '')')
         WRITE(LOU,'(''TITLE 6.5 0.8 SIZE 4 ''''|p| (GeV)'')')
C--
         WRITE(LOU,'(''('')')
         WRITE(LOU,'(''( AMST = '',1F7.3 )')  PTDATA(4,0,IE,0)
         WRITE(LOU,'(''( VTB2 = '',1F7.3 )')  PTDATA(6,0,IE,0)
         WRITE(LOU,'(''( AMSH = '',1E12.3)')  PTDATA(7,0,IE,0)
         WRITE(LOU,'(''( BETH = '',1F7.3 )')  PTDATA(8,0,IE,0)
         WRITE(LOU,'(''('')')
         WRITE(LOU,'(''SET ORDER X Y'')')
C--
         DO 2000 IA = 0, NA
            WRITE(LOU,'(''( E    = '',1F6.3)') PTDATA(3,0,IE,IA)
            WRITE(LOU,'(''( ALFS = '',1F6.3)') PTDATA(5,0,IE,IA)
            DO 200 IQ = 0, NQ
               Q    = PTDATA(1,IQ,IE,IA)
               DSG  = PTDATA(2,IQ,IE,IA)
               WRITE(LOU,*) Q, DSG
200         CONTINUE
            WRITE(LOU,'(''JOIN '',A8)') JOIN(IA)
2000     CONTINUE
20000 CONTINUE
C--
C  That's it.
C--
      STOP
      END
