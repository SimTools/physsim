C* (Update Record)
C*   96/07/23  K.Fujii		Added Charm and Tau modes.
C*				TAUON = 0 switches off tau mode.
C*
      SUBROUTINE INHGCN
 
      IMPLICIT    REAL*4 ( A-H, O-Z    )
C--
#include "../hg/include/HGCONS.inc"
#include "../sm/include/SMCONS.inc"
#include "../sm/include/SMPTAB.inc"
#include "../ss/include/SSCONS.inc"
C--
      REAL*8   AMA8, TNBET8, AMSQ8, AMSU8, AT8, AB8, XMU8,
     .         AMASS8(6), ALPHA8
      REAL*8   GV2PB8, ALF08, ALF8, S2W8, AMZ8, GMZ8, GMW8, AMH8, AMT8
      EXTERNAL BLKDAT
C--
      DATA AMSUSY / 1000. / 
C
C========< Entry Point >================================================
C
C--
C  Initialize Higgs masses and mixing angles.
C--
      GV2PB8 = xGV2PB
      ALF08  = 1/137.035989D0
      ALF8   = xALF
      S2W8   = xSIN2W
      AMZ8   = AMZ
      GMZ8   = GMZTOT
      GMW8   = GMWTOT
      AMH8   = AMH
      AMT8   = AMF(3,1,2)
C--
C  Local R*4 parameters.
C--
      AMSTU  = AMF(3,2,1)
      AMSS   = AMF(2,2,2)
      AMSC   = AMF(2,1,2)
      AMSB   = AMF(3,2,2)
      AMST   = AMF(3,1,2)
      AMSZ   = AMZ
      AMSW   = AMW
      COLOR  = 3
      TAUON  = 1
C>>>
C      TAUON  = 0
C      AMSB   = 5
C>>>
C--      
      AMA8   = xAMA
      TNBET8 = xTANB
C--
      AMSQ8  = AMSUSY
      AMSU8  = AMSUSY
      AT8    = AMSUSY
      AB8    = AMSUSY
      XMU8   = xAMU
C--
      CALL EWINIT(GV2PB8,ALF08,ALF8,S2W8,AMZ8,GMZ8,GMW8,AMH8,AMT8)
      CALL MHIGGS(AMA8, TNBET8, AMT8, AMSQ8, AMSU8, AT8, AB8, XMU8,
     .            AMASS8 , ALPHA8 )
C--
C  Set /HGCONS/.
C--
      AMA    = xAMA
      AMHL   = AMASS8(1)
      AMHH   = AMASS8(2)
      AMHM   = SQRT(AMA**2+AMW**2)
      BETA   = ATAN(TNBET8)
      SNALF  = SIN(ALPHA8)
      CSALF  = COS(ALPHA8)
      TNBET  = TNBET8
      SNBET  = SIN(BETA)
      CSBET  = COS(BETA)
      CSAMB  = CSALF*CSBET + SNALF*SNBET
      SNAMB  = SNALF*CSBET - CSALF*SNBET
C--
C  Print out SUSY parameters.
C--
      PRINT *, '             '
      PRINT *, ' ****** INHGCN uses the following SUSY parameters. ***'
      PRINT *, '             '
      PRINT *, '    AMSQ   = ', REAL(AMSQ8)
      PRINT *, '    AMSU   = ', REAL(AMSU8)
      PRINT *, '    AT     = ', REAL(AT8)
      PRINT *, '    AB     = ', REAL(AB8)
      PRINT *, '    XMU    = ', REAL(XMU8)
      PRINT *, '             '
      PRINT *, '    SNALF  = ', SNALF
      PRINT *, '    CSALF  = ', CSALF
      PRINT *, '    SNBET  = ', SNBET
      PRINT *, '    CSBET  = ', CSBET
      PRINT *, '    CSAMB  = ', CSAMB
      PRINT *, '    SNAMB  = ', SNAMB
C--
C  Initialize lighter Higgs width.
C--
      CALL GMH2FF(AMHL,AMSTU,AMSZ,GMH2TU)
      CALL GMH2FF(AMHL,AMSC,AMSZ,GMH2C)
      CALL GMH2FF(AMHL,AMSB,AMSZ,GMH2B)
      CALL GMH2FF(AMHL,AMST,AMSZ,GMH2T)
      CALL GMH2VV(AMHL,AMSW,AMSZ,GMH2W)
      CALL GMH2VV(AMHL,AMSZ,AMSZ,GMH2Z)
      C      = COLOR
      GMH2TU = TAUON*GMH2TU*(SNALF/CSBET)**2
      GMH2C  = C*GMH2C*(CSALF/SNBET)**2
      GMH2B  = C*GMH2B*(SNALF/CSBET)**2
      GMH2T  = C*GMH2T*(CSALF/SNBET)**2
      GMH2W  = 2*GMH2W*(SIN(ALPHA8-BETA))**2
      GMH2Z  = GMH2Z*(SIN(ALPHA8-BETA))**2
      GMHL   = GMH2TU + GMH2C + GMH2B + GMH2T + GMH2W + GMH2Z
      AMH    = AMHL
      GMHTOT = GMHL
C--
      PRINT *, '        '
      PRINT *, '    AMHL   = ', AMHL
      PRINT *, '    GMHL   = ', GMHL
      PRINT *, '    BR(tau)= ', GMH2TU/GMHL
      PRINT *, '    BR(c)  = ', GMH2C/GMHL
      PRINT *, '    BR(b)  = ', GMH2B/GMHL
      PRINT *, '    BR(t)  = ', GMH2T/GMHL
      PRINT *, '    BR(W)  = ', GMH2W/GMHL
      PRINT *, '    BR(Z)  = ', GMH2Z/GMHL
C--
C  Initialize heavier Higgs width.
C--
      CALL GMH2FF(AMHH,AMSTU,AMSZ,GMH2TU)
      CALL GMH2FF(AMHH,AMSC,AMSZ,GMH2C)
      CALL GMH2FF(AMHH,AMSB,AMSZ,GMH2B)
      CALL GMH2FF(AMHH,AMST,AMSZ,GMH2T)
      CALL GMH2VV(AMHH,AMSW,AMSZ,GMH2W)
      CALL GMH2VV(AMHH,AMSZ,AMSZ,GMH2Z)
      C      = COLOR
      GMH2TU = TAUON*GMH2TU*(CSALF/CSBET)**2
      GMH2C  = C*GMH2C*(SNALF/SNBET)**2
      GMH2B  = C*GMH2B*(CSALF/CSBET)**2
      GMH2T  = C*GMH2T*(SNALF/SNBET)**2
      GMH2W  = 2*GMH2W*(COS(ALPHA8-BETA))**2
      GMH2Z  = GMH2Z*(COS(ALPHA8-BETA))**2
      GMHH   = GMH2TU + GMH2C + GMH2B + GMH2T + GMH2W + GMH2Z
C--
      PRINT *, '        '
      PRINT *, '    AMHH   = ', AMHH
      PRINT *, '    GMHH   = ', GMHH
      PRINT *, '    BR(tau)= ', GMH2TU/GMHH
      PRINT *, '    BR(c)  = ', GMH2C/GMHH
      PRINT *, '    BR(b)  = ', GMH2B/GMHH
      PRINT *, '    BR(t)  = ', GMH2T/GMHH
      PRINT *, '    BR(W)  = ', GMH2W/GMHH
      PRINT *, '    BR(Z)  = ', GMH2Z/GMHH
C--
C  Initialize pseudo-scalar Higgs width.
C--
      CALL GMA2FF(AMA,AMSTU,AMSZ,GMH2TU)
      CALL GMA2FF(AMA,AMSC,AMSZ,GMH2C)
      CALL GMA2FF(AMA,AMSB,AMSZ,GMH2B)
      CALL GMA2FF(AMA,AMST,AMSZ,GMH2T)
      C      = COLOR
      GMH2TU = TAUON*GMH2TU*TNBET**2
      GMH2C  = C*GMH2C/TNBET**2
      GMH2B  = C*GMH2B*TNBET**2
      GMH2T  = C*GMH2T/TNBET**2
      GMA    = GMH2TU + GMH2C + GMH2B + GMH2T
C--
      PRINT *, '        '
      PRINT *, '    AMA    = ', AMA
      PRINT *, '    GMA    = ', GMA
      PRINT *, '    BR(tau)= ', GMH2TU/GMA
      PRINT *, '    BR(c)  = ', GMH2C/GMA
      PRINT *, '    BR(b)  = ', GMH2B/GMA
      PRINT *, '    BR(t)  = ', GMH2T/GMA
C--
C  Initialize charged Higgs width.
C--
      VTB  = VKM(3,3,2)
      VCB  = VKM(2,3,2)
      VCS  = VKM(2,2,2)
      CALL GMH2UD(AMHM,0.,AMSTU,AMSZ,TNBET,GMHTU)
      CALL GMH2UD(AMHM,AMST,AMSB,AMSZ,TNBET,GMHTB)
      CALL GMH2UD(AMHM,AMSC,AMSB,AMSZ,TNBET,GMHCB)
      CALL GMH2UD(AMHM,AMSC,AMSS,AMSZ,TNBET,GMHCS)
      C      = COLOR
      GMHTU  = TAUON*GMHTU
      GMHTB  = C*GMHTB*VTB**2
      GMHCB  = C*GMHCB*VCB**2
      GMHCS  = C*GMHCS*VCS**2
      GMHM   = GMHTU + GMHTB + GMHCB + GMHCS
C--
      PRINT *, '        '
      PRINT *, '    AMH-   = ', AMHM
      PRINT *, '    GMH-   = ', GMHM
      PRINT *, '    BR(tau)= ', GMHTU/GMHM
      PRINT *, '    BR(tb) = ', GMHTB/GMHM
      PRINT *, '    BR(cb) = ', GMHCB/GMHM
      PRINT *, '    BR(cs) = ', GMHCS/GMHM
C--
C  Modifiy top width by including t --> bH+ mode.
C--
      CALL GMT2HB(AMST,AMSB,AMHM,AMSW,TNBET,GMTHB)
      GMTBW      = GMF(3,1,2)
      GMTTOT     = GMTBW + GMTHB
      GMF(3,1,2) = GMTTOT
C--
      PRINT *, '        '
      PRINT *, '    AMST   = ', AMST
      PRINT *, '    GMTTOT = ', GMF(3,1,2)
      PRINT *, '    BR(bW) = ', GMTBW/GMTTOT
      PRINT *, '    BR(bH) = ', GMTHB/GMTTOT
      PRINT *, '        '
      PRINT *, ' *****************************************************'
C--
C  That's it.
C--
      RETURN
      END


