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
CINCLUDE DSGDP0
CINCLUDE GRQDSG
CINCLUDE GRQINT
CINCLUDE GRQADJ
CINCLUDE GRQPTP
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
C* Main program to draw the p_peak as a function of DE for
C* various alpha_s values.
C* (Update Record)
C*   93/04/20  K.Fujii       Full o(alpha_s) version.
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      PARAMETER    ( MXxPT = 200 )
      REAL     *8   PDATA(2,0:MXxPT)
      CHARACTER*8   JOIN(0:3)
C--
      DATA EPKSTD   / -2.5 /
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
C     AMST   = 150.D0
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
C  Initialize energy range.
C--
      NDE    = 10
      DEMN   = -1.D0
      DEMX   = +4.D0
      DDE    = (DEMX-DEMN)/NDE
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
C     CALL UALCPS(LOU,'TKSF.@@.@PTPK@.ALFS.TDR','RENEW','FB',30,IRT)
C--
C  Prepare TOPDRAW files.
C--
      WRITE(LOU,'(''( =========='')')
      WRITE(LOU,'(''(  p_peak '')')
      WRITE(LOU,'(''( =========='')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''NEWFRAME'')')
      WRITE(LOU,'(''SET FONT DUPLEX'')')
      WRITE(LOU,'(''SET INTENSITY 4'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET WINDOW X  2.5 12.5 Y 1.8  9.2'')')
      WRITE(LOU,'(''SET LIMIT X '',2E15.8)') DEMN, DEMX
      WRITE(LOU,'(''SET LIMIT Y 5.  25.'')')
      WRITE(LOU,'(''(SET SCALE Y LOG'')')
      WRITE(LOU,'(''SET LABELS SIZE 3.0'')')
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''TITLE 0.5 5.0 ANGLE 90 SIZE 4 ''''|p|0peak1'')')
      WRITE(LOU,'(''CASE                          ''''   X    X'')')
      WRITE(LOU,'(''TITLE 6.5 0.8 SIZE 4          ''''DE (GeV) '')')
      WRITE(LOU,'(''CASE                          ''''F        '')')
C--
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''( AMST = '',1F7.3 )')  AMST
      WRITE(LOU,'(''( VTB2 = '',1F7.3 )')  VTB2
      WRITE(LOU,'(''( AMSH = '',1E12.3)')  AMSH
      WRITE(LOU,'(''( BETH = '',1F7.3 )')  BETH
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET ORDER X DUMMY Y DUMMY'')')
C--
C  Loop over alpha_s.
C--
      MODE = 0
      DO 10000 IA = 0, NA
         ALPS  = ALPSMN + IA*DA
         PRINT *, ' alpha_s  = ', ALPS
C--
C  Hunt 1S peak.
C--
         EPK = EPKSTD
         CALL EPKHNT( 0, ALP,SN2W,AMSZ,AMSW,AMSB,
     .                AMST,ALPS,VTB2,AMSH,BETH, EPK,SGPK, NTRY )
         PRINT *, ' NTRY, EPK, SGPK =', NTRY, EPK, SGPK
C--
         WRITE(LOU,'(''( ALFS = '',1F6.3)') ALPS
         WRITE(LOU,'(''( EPK  = '',1F6.3)') EPK
C--
C  Loop over energy.
C--
         DO 1000 IDE = 0, NDE
            DE     = DEMN + IDE*DDE
            E      = EPK + DE
            RS     = 2*AMST +  E
            PRINT *, ' DE = ', DE, ' E = ', E
C--
C  |p| distribution.
C--
            SG    = 0
            DSMX  = 0
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
               CALL DSGDP0(MODE,ALPS,ALP,SN2W,AMSZ,AMSW,AMSB,
     .                     AMST, VTB2,  E,Q,DSG)
               MODE = 4
               SG   = SG + DSG*DQ*WT/3
C--
               PDATA(1,IQ) = Q
               PDATA(2,IQ) = DSG
C--
               IF ( DSG.GT.DSMX ) THEN
                  DSMX = DSG
                  IQMX = IQ
               ENDIF
100         CONTINUE
            MODE  = 3
            PRINT *, ' sqrt(s), sig_tt = ', RS, SG, ' pb'
C--
            PL    = PDATA(1,IQMX-1)
            PC    = PDATA(1,IQMX  )
            PH    = PDATA(1,IQMX+1)
            DSL   = PDATA(2,IQMX-1)
            DSC   = PDATA(2,IQMX  )
            DSH   = PDATA(2,IQMX+1)
C--
            A2    = (DSH+DSL)/2 - DSC
            A1    = (DSH-DSL)/2
            A0    = DSC
            XC    = -(A1/2/A2)
            PPK   = XC*(PC-PL) + PC
            DSGPK = A2*XC**2 + A1*XC + A0
C--
            WRITE(LOU,'(4E14.5)') DE, E, PPK, DSGPK
            WRITE(  *,'(4E14.5)') DE, E, PPK, DSGPK
1000     CONTINUE
         WRITE(LOU,'(''JOIN '',A8)') JOIN(IA)
         MODE  = 1
10000 CONTINUE
C--
C  That's it.
C--
      STOP
      END
