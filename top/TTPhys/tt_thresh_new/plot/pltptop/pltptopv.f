C*
C* Main program to draw the dsig/dp for various VTB2.
C* (Update Record)
C*   94/03/30  K.Fujii       Full o(alpha_s) version.
C*                           This version plots ptop for given DE's.
C*   94/08/03  K.Fujii       E1S is determined for each VTB2, and
C*                           is VTB2-dependent in this version.
C*
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      PARAMETER     ( MXxNP = 200, MXxNE = 20, MXxNX = 10 )
      REAL     *8   PTDATA(10,0:MXxNP,0:MXxNE,0:MXxNX)
C--
C  E1S is an approximate 1S peak position.
C--
C>>>
      REAL     *8   E1S
      PARAMETER     ( NE = 2 )
      REAL     *8   DEDATA(0:NE)
C     PARAMETER     ( NE = 0 )
C     REAL     *8   DEDATA(0:NE)
C>>>
C--
      CHARACTER*8   JOIN(0:4)
C>>>
      DATA E1S      / -2.5D0 /
      DATA DEDATA   /   0.D0, 1.D0, 2.D0 /
C     DATA DEDATA   /   2.D0 /
C>>>
C--
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT', ' ' /
      DATA LOU      / 21 /
C
C========< Entry Point >================================================
C
C--
C  Set parameters.
C--
      ALPS   = 0.12D0
      VTB2MN = 0.60D0
      VTB2MX = 1.40D0
      NV     = 4
      DV     = (VTB2MX-VTB2MN)/NV
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
C     CALL UALCPS(LOU,'TKSF.@@.@PTOP@.VTB2.IR.TDR','RENEW','FB',30,IRT)
C     CALL UALCPS(LOU,'TKSF.@@.@PTOP@.VTB2.BM.TDR','RENEW','FB',30,IRT)
C--
C  Loop over alpha_s.
C--
      MODE = 0
      DO 10000 IV = 0, NV
         VTB2  = VTB2MN + IV*DV
         PRINT *, ' VTB2 = ', VTB2
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
C--
               PTDATA(1,IQ,IE,IV) = Q
               PTDATA(2,IQ,IE,IV) = DSG
               PTDATA(3,IQ,IE,IV) = E
               PTDATA(4,IQ,IE,IV) = AMST
               PTDATA(5,IQ,IE,IV) = ALPS
               PTDATA(6,IQ,IE,IV) = VTB2
               PTDATA(7,IQ,IE,IV) = AMSH
               PTDATA(8,IQ,IE,IV) = BETH
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
         WRITE(LOU,'(''( ALFS = '',1F7.3 )')  PTDATA(5,0,IE,0)
         WRITE(LOU,'(''( AMSH = '',1E12.3)')  PTDATA(7,0,IE,0)
         WRITE(LOU,'(''( BETH = '',1F7.3 )')  PTDATA(8,0,IE,0)
         WRITE(LOU,'(''('')')
         WRITE(LOU,'(''SET ORDER X Y'')')
C--
         DO 2000 IV = 0, NV
            WRITE(LOU,'(''( E    = '',1F6.3)') PTDATA(3,0,IE,IV)
            WRITE(LOU,'(''( VTB2 = '',1F6.3)') PTDATA(6,0,IE,IV)
            DO 200 IQ = 0, NQ
               Q    = PTDATA(1,IQ,IE,IV)
               DSG  = PTDATA(2,IQ,IE,IV)
               WRITE(LOU,*) Q, DSG
200         CONTINUE
            WRITE(LOU,'(''JOIN '',A8)') JOIN(IV)
2000     CONTINUE
20000 CONTINUE
C--
C  That's it.
C--
      STOP
      END
