C*
C* Main program to draw the dsig/dp for various VTB2.
C* (Update Record)
C*   93/04/20  K.Fujii       Full o(alpha_s) version.
C*   94/08/03  K.Fujii       This version calculates E1S for
C*                           each VTB2 value.
C*
      IMPLICIT  REAL*8 ( A-H, O-Z )
      PARAMETER     ( MXxNP = 200, MXxNE = 50, MXxNX = 4 )
      REAL     *8   PTDATA(10,0:MXxNP,0:MXxNE,0:MXxNX)
      CHARACTER*8   JOIN(0:3)
C--
      DATA EPKSTD   / -2.5D0 /
      DATA VTBSTD   /  1.0D0 /
      DATA JOIN     / 'DOTDASH', 'SOLID', 'DASH', 'DOT' /
      DATA LOU      / 21 /
C
C========< Entry Point >================================================
C
C--
C  Set parameters.
C--
      ALPS   = 0.12D0
C>>>
      NV     = 2
      VTB2MN = 0.80D0
      VTB2MX = 1.20D0
C>>>
C     NV     = 1
C     VTB2MN = 0.60D0
C     VTB2MX = 1.40D0
C>>>
      DV     = (VTB2MX-VTB2MN)/NV
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
C  Initialize radiation function.
C--
      CALL ISRINT(0,2*AMST)
C--
C  Initialize energy range.
C--
      NDE    = 10
      DEMN   = -1.D0
      DEMX   = +4.D0
      DDE    = (DEMX-DEMN)/NDE
C>>>
C     DEMN   = +2.0D0
C     NDE    = 0
C>>>
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
C     CALL UALCPS(LOU,'TKSF.@@.@PTPK@.VTB2.IR.TDR','RENEW','FB',30,IRT)
C      CALL UALCPS(LOU,'TKSF.@@.@PTPK@.VTB2.BM.TDR','RENEW','FB',30,IRT)
C--
C  Loop over |V_tb|^2.
C--
      MODE = 0
      DO 10000 IV = 0, NV
         VTB2  = VTB2MN + IV*DV
         PRINT *, ' |V_tb|^2 = ', VTB2
C--
C  Hunt 1S peak.
C--
         EPK  = EPKSTD
         CALL EPKHNT( 0, ALP,SN2W,AMSZ,AMSW,AMSB,
     .                AMST,ALPS,VTB2,AMSH,BETH, EPK,SGPK, NTRY )
         PRINT *, ' NTRY, EPK, SGPK =', NTRY, EPK, SGPK
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
               PTDATA(1,IQ,IDE,IV) = Q
               PTDATA(2,IQ,IDE,IV) = DSG
               PTDATA(3,IQ,IDE,IV) = DE
               PTDATA(4,IQ,IDE,IV) = AMST
               PTDATA(5,IQ,IDE,IV) = ALPS
               PTDATA(6,IQ,IDE,IV) = VTB2
               PTDATA(7,IQ,IDE,IV) = AMSH
               PTDATA(8,IQ,IDE,IV) = BETH
               PTDATA(9,IQ,IDE,IV) = EPK
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
      WRITE(LOU,'(''( AMST = '',1F7.3 )')  PTDATA(4,0,0,0)
      WRITE(LOU,'(''( ALFS = '',1F6.3 )')  PTDATA(5,0,0,0)
      WRITE(LOU,'(''( AMSH = '',1E12.3)')  PTDATA(7,0,0,0)
      WRITE(LOU,'(''( BETH = '',1F7.3 )')  PTDATA(8,0,0,0)
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''SET ORDER X DUMMY Y DUMMY DUMMY'')')
C--
      DO 20000 IV = 0, NV
         WRITE(LOU,'(''( VTB2 = '',1F6.3)') PTDATA(6,0,0,IV)
         WRITE(LOU,'(''( EPK  = '',1F6.3)') PTDATA(9,0,0,IV)
         VTB2 = PTDATA(6,0,0,IV)
         DO 2000 IDE = 0, NDE
            DE   = PTDATA(3,0,IDE,0)
            E    = PTDATA(9,0,0,IV) + DE
            DSMX = 0
            DO 200 IQ = 0, NQ
               Q    = PTDATA(1,IQ,IDE,IV)
               DSG  = PTDATA(2,IQ,IDE,IV)
               IF ( DSG.GT.DSMX ) THEN
                  DSMX = DSG
                  IQMX = IQ
               ENDIF
200         CONTINUE
C--
            DSL   = PTDATA(2,IQMX-1,IDE,IV)
            DSC   = PTDATA(2,IQMX  ,IDE,IV)
            DSH   = PTDATA(2,IQMX+1,IDE,IV)
            PL    = PTDATA(1,IQMX-1,IDE,IV)
            PC    = PTDATA(1,IQMX  ,IDE,IV)
            PH    = PTDATA(1,IQMX+1,IDE,IV)
C--
            A2    = (DSH+DSL)/2 - DSC
            A1    = (DSH-DSL)/2
            A0    = DSC
            XC    = -(A1/2/A2)
            PPK   = XC*(PC-PL) + PC
            DSGPK = A2*XC**2 + A1*XC + A0
C--
            WRITE(LOU,'(5E14.5)') DE, E, PPK, DSGPK, VTB2
            WRITE(*,'(5E14.5)') DE, E, PPK, DSGPK, VTB2
2000     CONTINUE
         WRITE(LOU,'(''JOIN '',A8)') JOIN(IV)
20000 CONTINUE
C--
C  That's it.
C--
      STOP
      END

