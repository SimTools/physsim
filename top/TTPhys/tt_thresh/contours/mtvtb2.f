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
CINCLUDE FHIGGS
CINCLUDE POTHIG
CINCLUDE POTQCD
C>>>
CINCLUDE GTCHI2V
CINCLUDE UCONT22
CINCLUDE USOLVE
CINCLUDE RSHDIS
CINCLUDE UDSRCH
C*
C* (Update Record)
C*   91/09/02  K.Fujii       Original version for Finland conference.
C*                           This version draws m_t-!Vtb!^2 contours.
C*
      IMPLICIT  REAL*8 ( A-H, O-Y )
C     PARAMETER   ( MXxAM = 20, MXxLM = 10, MXxCNT = 3 )
      PARAMETER   ( MXxAM = 40, MXxLM = 14, MXxCNT = 3 )
      REAL*8       RDATA(3,0:MXxAM,0:MXxLM), CHICNT(MXxCNT)
C--
      PARAMETER   ( MXxRS = 11 )
      REAL*8       EXPCND(4,MXxRS)
C--
C170      DATA AMTMN / 169.8D0 / AMTMX / 170.2D0 /
C170      DATA ALFMN / 0.510D0 / ALFMX / 1.490D0 /
      DATA AMTMN / 174.6D0 / AMTMX / 175.4D0 /
      DATA ALFMN / 0.300D0 / ALFMX / 1.700D0 /
      DATA (EXPCND(1,I),I=1,MXxRS) /  -7.D0, -6.D0, -5.D0,
     .                                -4.D0, -3.D0, -2.D0,
     .                                -1.D0,  0.D0,  1.D0,
     .                                 2.D0,  3.D0/
      DATA (EXPCND(2,I),I=1,MXxRS) / MXxRS*1.000D+3  /
      DATA (EXPCND(3,I),I=1,MXxRS) / MXxRS*0.260D+0  /
      DATA (EXPCND(4,I),I=1,MXxRS) / MXxRS*0.852D-2  /
C--
      DATA LOU /20/
C
C========< Entry Point >================================================
C
C--
C  Allocate output data set.
C--
C     CALL UALCPS(LOU,'TKSF.@.@MTVTB2@.CONT.TDR','RENEW','FB',30,IRT)
C--
C  Nominal m_t and lambda.
C--
C     AMT   = 150
C     AMT   = 170
      AMT   = 175
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
         VTB2 = ALMN + DLM*ILM
         MODE = 0
         DO 10 IAM = 0, NAM
            AMT = AMN + DAM*IAM
            CALL GTCHI2(MODE,AMT,ALFS,VTB2,AMH,BTH,MXxRS,EXPCND,CHI2)
            RDATA(1,IAM,ILM) = VTB2
            RDATA(2,IAM,ILM) = AMT
            RDATA(3,IAM,ILM) = CHI2
C>>>
            WRITE(LOU,'(''(  '',3F15.6)')  AMT, VTB2, CHI2
C>>>
C           PRINT *, ' AMT, VTB2, CHI2 = ', AMT, VTB2, CHI2
C>>>
            IF ( CHI2.LT.CHI2MN ) THEN
               XMN    = AMT
               YMN    = VTB2
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
C     PARAMETER    ( MXxRS = 10, MXxAM = 20, MXxAL = 10 )
      PARAMETER    ( MXxRS = 10, MXxAM = 40, MXxAL = 14 )
C>>>
      REAL*8       SGDAT(4,0:MXxRS,0:MXxAM,0:MXxAL)
      REAL*8       SGOBS(2,MXxRS+1)
C>>>  For MT-ALFS
C     DATA ISEED   / 31414927 /
C     DATA ISEED   / 1        /
C>>>  For MT-VTB2
      DATA ISEED   / 3        /
C>>>  For MH-BTH2
C     DATA ISEED   / 31414927 /
C>>>
      DATA ZERO    / 1.D-4  /
C     DATA LIU     / 8 /
      DATA LIU     / 5 /
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
            PRINT *, ' RS = ', RS, ' SGTH = ', SGD,
     .               ' SGEXP = ', SGOBS(1,IRS), '+/-', SGOBS(2,IRS)
            WRITE(20,'(''('',4E15.5)') RS, SGOBS(1,IRS), SGOBS(2,IRS),
     .                                 ENOBS
C>>>
            CHI2  = CHI2 + ((SGD-SGOBS(1,IRS))/SGOBS(2,IRS))**2
            MODEP = 3
10       CONTINUE
C--
C  Read in cross section data.
C--
         PRINT *, 'Reading cross section data'
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
      IAL   = (VTB2-ALMN+ZERO)/DAL
      IAM   = (AMT -AMMN+ZERO)/DAM
      CHI2  = 0
      DO 30 IRS = 1, NRS
         RS    = EXPCND(1,IRS)
         JRS   = (RS-RSMN+ZERO)/DRS
         ALFSD = SGDAT(1,JRS,IAM,IAL)
         AMTD  = SGDAT(2,JRS,IAM,IAL)
         RSD   = SGDAT(3,JRS,IAM,IAL)
         SGD   = SGDAT(4,JRS,IAM,IAL)
         IF ( ABS(RS-RSD).GT.1.D-4 .OR.
     .        ABS(AMT-AMTD).GT.1.D-4 .OR.
     .        ABS(VTB2-ALFSD).GT.1.D-4 ) THEN
            PRINT *, ' ALFS = ', VTB2, ' ALFSD = ', ALFSD
            PRINT *, ' AMT  = ', AMT , ' AMTD  = ', AMTD
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
