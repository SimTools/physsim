CC**********************************************************************
C*
C*====================
C* Subroutine SPPARM
C*====================
C*
C* (Function)
C*   Check if inputs from FT05 are consistent with those used in the
C*   integration step, when NGETBS = 1.
C*   When NGETBS = 0, write the inputs from FT05 to BASES:Inputs.
C* (Update Record)
C*   95/03/21  K.Fujii		Original version for TTBARZ.
C*   95/06/07  K.Fujii		Added ISRB switch.
C*
CC**********************************************************************

      SUBROUTINE SPPARM

#include "USRPRM.inc"
#include "USMPRM.inc"
#include "BS51PR.inc"
C--
      PARAMETER      ( MX_WRK = 1000 )
      COMMON /SSWORK/ IDATA
      INTEGER*4       IDATA(MX_WRK)
      REAL   *4       RDATA(MX_WRK)
      EQUIVALENCE    (IDATA,RDATA)
      CHARACTER*16    GNNAME
      EQUIVALENCE    (IDATA(3),GNNAME)
C--      
      CHARACTER*44 BKNAME / 'BASES:Inputs' /
C
C =====< Entry Point >==================================================
C
C--
C  Branch on NGETBS.
C--
      IF ( NGETBS.EQ.0 ) THEN
C--
C  Need to do BASES. TBCRTE 'BASES:Inputs'.
C--
         IF ( NPUTBS.EQ.1 ) THEN
            PRINT *, ' '
            PRINT *, 'Write the inputs from FT05 to ', BKNAME
            PRINT *, ' '
            CALL TBCRTE(2,BKNAME, 0, 0, IRT)
            IF ( IRT.LT.0 ) THEN
               PRINT *, ' >>>SPPARM: TBCRTE IRT=', IRT, ' for ', BKNAME
               PRINT *, '    STOP!'
               STOP
            ENDIF
C--
C  TBPUT element #1 of 'BASES:Inputs'.
C--
            CALL UIDATE(IDATE)
            ITIME    = TIME()
            IDATA(1) = IDATE
            IDATA(2) = ITIME
            GNNAME   = 'ttz'
            IDATA(7) = 0205
C--
            IE = 1
            NW = 7
            CALL TBPUT(2,BKNAME,IE,NW,IDATA,IRT)
            IF ( IRT.LT.0 ) THEN
               PRINT *, ' >>>SPPARM: TBPUT IRT=', IRT, ' for element #',
     .                  IE, ' of ', BKNAME               
               PRINT *, '    STOP!'
               STOP
            ENDIF
C--
C  TBPUT element #2 of 'BASES:Inputs'.
C--
            RDATA(1) = ALFI
            RDATA(2) = ALFS
            RDATA(3) = AMSW
            RDATA(4) = AMSZ
            RDATA(5) = AMSH
            RDATA(6) = AMST
C--
            IE = 2
            NW = 6
            CALL TBPUT(2,BKNAME,IE,NW,RDATA,IRT)
            IF ( IRT.LT.0 ) THEN
               PRINT *, ' >>>SPPARM: TBPUT IRT=', IRT, ' for element #',
     .                  IE, ' of ', BKNAME               
               PRINT *, '    STOP!'
               STOP
            ENDIF
C--
C  TBPUT element #4 of 'BASES:Inputs'.
C--
            RDATA(1) = SQRTS
            RDATA(2) = POLEBM
            RDATA(3) = SGMEBM
            RDATA(4) = 0
            IDATA(5) = ISRB
C--
            IE = 4
            NW = 5
            CALL TBPUT(2,BKNAME,IE,NW,IDATA,IRT)
            IF ( IRT.LT.0 ) THEN
               PRINT *, ' >>>SPPARM: TBPUT IRT=', IRT, ' for element #',
     .                  IE, ' of ', BKNAME               
               PRINT *, '    STOP!'
               STOP
            ENDIF
         ENDIF
      ELSE IF ( NGETBS.EQ.1 ) THEN
C--
C  Input parameters should be replaced by those in BASES:Inputs.
C--
         PRINT *, ' '
         PRINT *, 'SPPARM: going to replace inputs from FT05 by ',
     .            'those used in the integration step.'
         PRINT *, ' '     
C--
C  TBGET element #2 of 'BASES:Inputs'.
C--
         IE = 2
         CALL TBGET(2,BKNAME,IE,NW,IDATA,IRT)
         IF ( IRT.LT.0 ) THEN
            PRINT *, ' >>>SPPARM: TBGET IRT=', IRT, ' for element #',
     .               IE, ' of ', BKNAME               
            PRINT *, '    STOP!'
            STOP
         ENDIF
C--
         ALFI = RDATA(1)
         ALFS = RDATA(2)
         AMSW = RDATA(3)
         AMSZ = RDATA(4)
         AMSH = RDATA(5)
         AMST = RDATA(6)
C--
C  TBGET element #4 of 'BASES:Inputs'.
C--
         IE = 4
         CALL TBGET(2,BKNAME,IE,NW,IDATA,IRT)
         IF ( IRT.LT.0 ) THEN
            PRINT *, ' >>>SPPARM: TBGET IRT=', IRT, ' for element #',
     .               IE, ' of ', BKNAME               
            PRINT *, '    STOP!'
            STOP
         ENDIF
C--
         SQRTS  = RDATA(1)
         POLEBM = RDATA(2)
         SGMEBM = RDATA(3)
         ISRB   = IDATA(5)
C--
C  Print modified user parameters.
C--
         PRINT *,' '
         PRINT *,'SPPARM Modified User Parameters'
         PRINT *,'-------------------------------'
         PRINT *,' '
         PRINT *,' ALFI=',ALFI,' ALFS=',ALFS
         PRINT *,' AMSW=',AMSW,' AMSZ=',AMSZ,' AMSH=',AMSH
         PRINT *,' AMST=',AMST
         PRINT *,' '
         PRINT *,' ECM =',SQRTS
         PRINT *,' POLE=',POLEBM
         PRINT *,' SGEB=',SGMEBM
         PRINT *,' ISRB=',ISRB
         PRINT *,' '
      ELSE 
C--
C  Input parameters will not be touched.
C--
         PRINT *, 'Input parameters will not be touched by SPPARM.'
         PRINT *, 'Users should make sure that inputs in the BASES'
         PRINT *, 'section are consistent with those used in the '
         PRINT *, 'integration step.'
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
