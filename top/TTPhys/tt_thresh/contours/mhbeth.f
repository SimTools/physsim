CINCLUDE SGTTEF
CINCLUDE SGTTHR2
CINCLUDE RKUTTA2
CINCLUDE EI
CINCLUDE FX
CINCLUDE FH
CINCLUDE F1
CINCLUDE GTGAMT
C>>>
C For o(alpha_s) correction
C*INCLUDE SETPRM
CINCLUDE SETPRMC
CINCLUDE TBWCORR
C>>>
C*INCLUDE POTRCH
C*INCLUDE POTCPL
C*INCLUDE POTMTN
CINCLUDE POTHIG
CINCLUDE FHIGGS
CINCLUDE POTQCD
C>>>
CINCLUDE GTCHI2H
CINCLUDE UCONT22
CINCLUDE USOLVE
CINCLUDE RSHDIS
CINCLUDE UDSRCH
C*
C* (Update Record)
C*   91/09/02  K.Fujii       Original version for Finland conference.
C*                           This version draws m_H-beta_H^2 contours.
C*
      IMPLICIT  REAL*8 ( A-H, O-Y )
      PARAMETER   ( MXxAM = 10, MXxLM = 6, MXxCNT = 3 )
      REAL*8       RDATA(3,0:MXxAM,0:MXxLM), CHICNT(MXxCNT)
C--
      PARAMETER   ( MXxRS = 11 )
      REAL*8       EXPCND(4,MXxRS)
C--
      DATA AMTMN /  50.0D0 / AMTMX / 300.0D0 /
      DATA ALFMN / 1.00D-3 / ALFMX / 4.000D0 /
      DATA (EXPCND(1,I),I=1,MXxRS) /  -7.D0, -6.D0, -5.D0,
     .                                -4.D0, -3.D0, -2.D0,
     .                                -1.D0,  0.D0,  1.D0,
     .                                 2.D0,  3.D0/
C>>>
      DATA (EXPCND(2,I),I=1,MXxRS) / MXxRS*1.000D+3  /
C     DATA (EXPCND(2,I),I=1,MXxRS) / MXxRS*2.000D+3  /
      DATA (EXPCND(3,I),I=1,MXxRS) / MXxRS*0.260D+0  /
      DATA (EXPCND(4,I),I=1,MXxRS) / MXxRS*0.852D-2  /
C>>>
C--
      DATA LOU /20/
C
C========< Entry Point >================================================
C
C--
C  Allocate output data set.
C--
C     CALL UALCPS(LOU,'TKSF.@.@MHBETH@.CONT.TDR','RENEW','FB',30,IRT)
C--
C  Nominal m_t and lambda.
C--
C     AMT   = 150
      AMT   = 170
      ALFS  = 0.12D0
C--
C  Initialize sqrt(s).
C--
      DO 1 IRS = 1, MXxRS
         EXPCND(1,IRS) = EXPCND(1,IRS) + 2*AMT
1     CONTINUE
C--
C  Initialization.
C--
      AMN   = AMTMN
      AMX   = AMTMX
      NAM   = MXxAM
      DAM   = (AMX-AMN)/NAM
      ALMN  = ALFMN
      ALMX  = ALFMX
      NLM   = MXxLM
      DLM   = (ALMX-ALMN)/NLM
      VTB2  = 1
      AMH   = 1.D10
      BTH   = 1
C>>>
      PRINT *, ' First create MC data.'
C>>>
      CALL GTCHI2(-1,AMT,ALFS,VTB2,AMH,BTH,MXxRS,EXPCND,CHI2)
C>>>
      PRINT *, '    CHI2ST = ', CHI2
C--
C  Start scanning m_t-lambda mesh.
C--
C>>>
      PRINT *, ' Then create mesh data.'
C>>>
      CHI2MN = 1.D20
      DO 100 ILM = 0, NLM
         BTH2 = ALMN + DLM*ILM
         BTH  = SQRT(BTH2)
         MODE = 0
         DO 10 IAM = 0, NAM
            AMH = AMN + DAM*IAM
            CALL GTCHI2(MODE,AMT,ALFS,VTB2,AMH,BTH,MXxRS,EXPCND,CHI2)
            RDATA(1,IAM,ILM) = BTH2
            RDATA(2,IAM,ILM) = AMH
            RDATA(3,IAM,ILM) = CHI2
C>>>
            WRITE(LOU,'(''(  '',3F15.6)')  AMH, BTH2, CHI2
C>>>
C           PRINT *, ' AMH, BTH2, CHI2 = ', AMH, BTH2, CHI2
C>>>
            IF ( CHI2.LT.CHI2MN ) THEN
               XMN    = AMH
               YMN    = BTH2
               CHI2MN = CHI2
            ENDIF
            MODE = 1
10       CONTINUE
100   CONTINUE
C>>>
      PRINT *, ' XMN, YMN, CHI2MN = ', XMN, YMN, CHI2MN
      WRITE(LOU,'(''('')')
      WRITE(LOU,'(''( Chi^2 minimum = '',E15.5)') CHI2MN
      WRITE(LOU,'(''(  at X = '',E15.5,'' Y = '',E15.5)') XMN, YMN
      WRITE(LOU,'(''('')')
C--
C  Make countour.
C--
C>>>
      PRINT *, ' Now satrt drawing countour.'
C>>>
      CHI1S = CHI2MN + 1.00
      CHI68 = CHI2MN + 2.28
      CHI90 = CHI2MN + 4.61
      CHICNT(1) = CHI1S
      CHICNT(2) = CHI68
      CHICNT(3) = CHI90
C--
CX    CALL UCONT2(LOU,MXxCNT,CHICNT,MXxLM,MXxAM,RDATA)
      CALL UCONT2(LOU,MXxCNT,CHICNT,MXxAM,MXxLM,RDATA)
C--
C  That's it.
C--
      STOP
      END

C--
C     MODE = -1 : first call.
C          =  0 : if alpha_s changed.
C          =  1 : if alpha_s unchanged.
C--
 
      SUBROUTINE GTCHI2(MODE,AMT,ALFS,VTB2,AMH,BTH,NRS,EXPCND,CHI2)
 
      IMPLICIT     REAL*8 ( A-H, O-Y )
      INTEGER*4    MODE, NRS
      REAL*8       AMT, ALFS, VTB2, AMH, BTH, EXPCND(4,NRS), CHI2
C>>>
      PARAMETER    ( MXxRS = 10, MXxAM = 10, MXxAL = 6 )
C>>>
      REAL*8       SGDAT(4,0:MXxRS,0:MXxAM,0:MXxAL)
      REAL*8       SGOBS(2,MXxRS)
C>>>  For MT-ALFS
C     DATA ISEED   / 3        /
C>>>  For MT-VTB2
C     DATA ISEED   / 3        /
C>>>  For MH-BTH2
      DATA ISEED   / 3        /
C     DATA ISEED   / 31415927  /
C>>>
      DATA ZERO    / 1.D-4  /
      DATA LIU     / 8 /
C
C========< Entry Point >================================================
C
C--
C  First calculate nominal cross sections.
C     EXPCND(1,*) = sqrt(s)
C           (2,*) = integrated luminosity
C           (3,*) = acceptance
C           (4,*) = background cross section * efficiency
C--
      IF ( MODE.EQ.-1 ) THEN
         MODEP = 0
C>>>
         PRINT *, 'Create MC data with the following parameters.'
         PRINT *, '  AMT  = ', AMT
         PRINT *, '  ALFS = ', ALFS
         PRINT *, '  VTB2 = ', VTB2
         PRINT *, '  AMH  = ', AMH
         PRINT *, '  BTH  = ', BTH
C>>>
         CHI2 = 0
         DO 10 IRS = 1, NRS
            RS    = EXPCND(1,IRS)
            CALL SGTTEF(MODEP,AMT,ALFS,VTB2,AMH,BTH,RS,SGD)
            SGD   = SGD
            SLM   = EXPCND(2,IRS)
            ACC   = EXPCND(3,IRS)
            SGBG  = EXPCND(4,IRS)
            ENSG  = SLM*ACC*SGD
            ENBG  = SLM*SGBG
            ENOBS = ENSG + ENBG
            ENERR = SQRT(ENOBS)
            ENOBS = ENOBS + ENERR*RANN(ISEED)
            ENERR = SQRT(ENOBS)
            ENOBS = ENOBS - ENBG
            SGOBS(1,IRS) = ENOBS/SLM/ACC
            SGOBS(2,IRS) = ENERR/SLM/ACC
C>>>
C           PRINT *, ' RS = ', RS, ' SGTH = ', SGD,
C    .               ' SGEXP = ', SGOBS(1,IRS), '+/-', SGOBS(2,IRS)
            WRITE(20,'(''('',4E15.5)') RS, SGOBS(1,IRS), SGOBS(2,IRS),
     .                                 ENOBS
C>>>
            CHI2  = CHI2 + ((SGD-SGOBS(1,IRS))/SGOBS(2,IRS))**2
            MODEP = 3
10       CONTINUE
C--
C  Read in cross section data.
C--
         DO 2000 IAL = 0, MXxAL
            DO 200 IAM = 0, MXxAM
               DO 20 IRS = 0, MXxRS
                  READ(LIU,*) (SGDAT(I,IRS,IAM,IAL),I=1,4)
C                 PRINT *,    (SGDAT(I,IRS,IAM,IAL),I=1,4)
20             CONTINUE
200         CONTINUE
2000     CONTINUE
C--
         ALMN   = SGDAT(1,0,0,0)
         ALMX   = SGDAT(1,MXxRS,MXxAM,MXxAL)
         DAL    = (ALMX-ALMN)/MXxAL
         PRINT *, ' ALMN, ALMX, DAL = ', ALMN, ALMX, DAL
C--
         AMMN   = SGDAT(2,0,0,0)
         AMMX   = SGDAT(2,MXxRS,MXxAM,MXxAL)
         DAM    = (AMMX-AMMN)/MXxAM
         PRINT *, ' AMMN, AMMX, DAM = ', AMMN, AMMX, DAM
C--
         RSMN   = SGDAT(3,0,0,0)
         RSMX   = SGDAT(3,MXxRS,MXxAM,MXxAL)
         DRS    = (RSMX-RSMN)/MXxRS
         PRINT *, ' RSMN, RSMX, DRS = ', RSMN, RSMX, DRS
C--
         RETURN
      ENDIF
C--
C  Loop over energyies and accumulate chi**2.
C--
      BTH2  = BTH**2
      IAL   = (BTH2-ALMN+ZERO)/DAL
      IAM   = (AMH -AMMN+ZERO)/DAM
      CHI2  = 0
      DO 30 IRS = 1, NRS
         RS    = EXPCND(1,IRS)
         JRS   = (RS-RSMN+ZERO)/DRS
         ALFSD = SGDAT(1,JRS,IAM,IAL)
         AMTD  = SGDAT(2,JRS,IAM,IAL)
         RSD   = SGDAT(3,JRS,IAM,IAL)
         SGD   = SGDAT(4,JRS,IAM,IAL)
         IF ( ABS(RS-RSD).GT.1.D-4 .OR.
     .        ABS(AMH-AMTD).GT.1.D-4 .OR.
     .        ABS(BTH2-ALFSD).GT.1.D-4 ) THEN
            PRINT *, ' BTH2 = ', BTH2, ' BTH2D = ', ALFSD
            PRINT *, ' AMH  = ', AMH , ' AMHD  = ', AMTD
            PRINT *, ' RS   = ', RS,   ' RSD   = ', RSD
            STOP
         ENDIF
         CHI2  = CHI2 + ((SGD-SGOBS(1,IRS))/SGOBS(2,IRS))**2
30    CONTINUE
C--
C  That's it.
C--
      RETURN
      END
