C*
C* Main program to draw the dsig/dp for various AMST.
C* (Update Record)
C*   93/04/20  K.Fujii       Full o(alpha_s) version.
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      PARAMETER     ( MXxNP = 200, MXxNE = 20, MXxNX = 10 )
      REAL     *8   PTDATA(10,0:MXxNP,0:MXxNE,0:MXxNX)
C--
C  The first entry is a point below 1S peak.
C--
      PARAMETER     ( NE = 0 )
C     PARAMETER     ( NE = 2 )
      REAL     *8   EDATA(0:NE)
C--
      CHARACTER*8   JOIN(0:3)
C>>>
C     DATA EDATA    / -2.5D0, 0.D0, 1.D0 /
      DATA EDATA    / -2.5D0 /
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
C     AMSTMN = 149.D0
C     AMSTMX = 151.D0
      AMSTMN = 169.D0
      AMSTMX = 171.D0
C
      NM     = 2
      DM     = (AMSTMX-AMSTMN)/NM
      ALPS   = 0.12D0
      VTB2   = 1.D0
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
      CALL ISRINT(0,AMSTMN+AMSTMX)
C--
C  Delta_E for IE = 0.
C--
      DE     = +1.D0
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
C     CALL UALCPS(LOU,'TKSF.@@.@PTOP@.AMST.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@PTOP@.AMST.BM.TDR','RENEW','FB',30,IRT)
C--
C  Loop over alpha_s.
C--
      MODE  = 0
      MODEP = 0
      DO 10000 IM = 0, NM
         AMST  = AMSTMN + IM*DM
         PRINT *, ' m_t  = ', AMST
C--
C  Loop over energy.
C--
         DO 1000 IE = 0, NE
            E      = EDATA(IE)
            IF ( IE.EQ.0 ) THEN
               EPK = E
               CALL EPKHNT(MODEP, ALP,SN2W,AMSZ,AMSW,AMSB,
     .                    AMST,ALPS,VTB2,AMSH,BETH, EPK,SGPK, NTRY )
               PRINT *, ' NTRY, EPK, SGPK =', NTRY, EPK, SGPK
               E   = EPK + DE
            ENDIF
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
               PTDATA(1,IQ,IE,IM) = Q
               PTDATA(2,IQ,IE,IM) = DSG
               PTDATA(3,IQ,IE,IM) = E
               PTDATA(4,IQ,IE,IM) = AMST
               PTDATA(5,IQ,IE,IM) = ALPS
               PTDATA(6,IQ,IE,IM) = VTB2
               PTDATA(7,IQ,IE,IM) = AMSH
               PTDATA(8,IQ,IE,IM) = BETH
100         CONTINUE
            MODE  = 3
            PRINT *, ' sqrt(s), sig_tt = ', RS, SG, ' pb'
1000     CONTINUE
         MODE  = 2
         MODEP = 1
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
C--
         WRITE(LOU,'(''( =================='')')
         IF ( IE.EQ.0 ) THEN
            WRITE(LOU,'(''(  E = E1S'')')
         ELSE
            WRITE(LOU,'(''(  E = '',1F8.4)') E
         ENDIF
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
            WRITE(LOU,'(''TITLE 5.0 8.7 SIZE 3 ''''+'',1F4.1,''GeV'')')
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
         WRITE(LOU,'(''( ALFS = '',1F7.3 )')  PTDATA(5,0,IE,0)
         WRITE(LOU,'(''( VTB2 = '',1F7.3 )')  PTDATA(6,0,IE,0)
         WRITE(LOU,'(''( AMSH = '',1E12.3)')  PTDATA(7,0,IE,0)
         WRITE(LOU,'(''( BETH = '',1F7.3 )')  PTDATA(8,0,IE,0)
         WRITE(LOU,'(''('')')
         WRITE(LOU,'(''SET ORDER X Y'')')
C--
         DO 2000 IM = 0, NM
            WRITE(LOU,'(''( E    = '',1F6.3)') PTDATA(3,0,IE,IM)
            WRITE(LOU,'(''( AMST = '',1F6.1)') PTDATA(4,0,IE,IM)
            DO 200 IQ = 0, NQ
               Q    = PTDATA(1,IQ,IE,IM)
               DSG  = PTDATA(2,IQ,IE,IM)
               WRITE(LOU,*) Q, DSG
200         CONTINUE
            WRITE(LOU,'(''JOIN '',A8)') JOIN(IM)
2000     CONTINUE
20000 CONTINUE
C--
C  That's it.
C--
      STOP
      END
