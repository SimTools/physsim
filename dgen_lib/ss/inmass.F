      SUBROUTINE INMASS(ALF,ALFS,WM,ZM,AM0,AMU,AM1,AM2,AM3,TANB,AMA,IS,
     &           xSFM, SFM3,ASFT,FM3,
     &           SWM,SZM,SFM,SHM,ALPH,SFM3M,SPHI)
C
C     This program calculate MSUGRA mass spectrum from 
C     universal paramters
C     universal scalar mass AM0 
C     higgsino mass         AMU
C     gaugino mass          AM2
C     tanbeta               TANB
C     pseudoscaler mass     AMA (to calculate higgs sector)
C  
C     it also calculates neutralino mixing angle, 
C     and sfermion  mixing angle 
C
C     Notice if IS=0, program calcualtes weak scale masses,  
C     by simplified GUT relation  while if IS=2, it 
C     uses sfermion mass parameter and ino mass 
C     parameters (sfermion masses 
C     and AM1 AM2 AM3) taken from input.  
C     If IS=1, only third generation mass is the one taken 
C     from input
      IMPLICIT    REAL*8 ( A-H, O-Z )
#include "../ss/include/SSMIX.inc"
      REAL*8   OL(2,2), OR(2,2)
      REAL*8   FM(3), FM3(3), xSFM(7),SFM3(7), ASFT(3)
      REAL*8   SWM(2), SZM(4),SFM(7),SHM(4),SFM3M(7),SPHI(3)
      REAL*8   ALF, ALFS, ALPH, WM, ZM, S2W, M1,M2,M3
     .         AM0, AMU, AM1, AM2, TANB, AMA
      COMPLEX*16 UN(4,4)
C
 
C Initialize  SSINPT.inc
C--
      IF(IS.EQ.2) THEN 
       DO I=1,7
       SFM(I)=xSFM(I)
       ENDDO
       M1=AM1
       M2=AM2
       M3=AM3
      ENDIF
      
C       DO I=1,7
C       SFM3(I)=xSFM3(I)
C       ENDDO
C       DO I=1,3
C       ASFT(I)=xASFT(I)
C       ENDDO

C--
C  Initialize particle table and coupling constants.
C--
      S2W   = ( 1 - WM/ZM )*( 1 + WM/ZM )
      IF(IS.NE.2) THEN 
C GUT RELATION TO GAUGINO MASS 
       M2=AM2
       M1=S2W/(1.D0-S2W)*M2/0.6D0
       AM1 = M1
       M3=ALFS*S2W/ALF*M2
      ENDIF
      FM(1) = 0
      FM(2) = 0
      FM(3) = 0
cccccccccc test at 11/16/98 cccccccccc
ccc      M1=87.83
ccc      M1 = 87.83585
      write(6,*) '========== INMASS =========='
      write(6,*) 'IS=',IS
      write(6,*) 'M0=',AM0
      write(6,*) 'M1=',M1
      write(6,*) 'M2=',M2
      write(6,*) 'M3=',M3
      write(6,*) 'SFM3(1) = ',SFM3(1)
      write(6,*) 'SFM3(2) = ',SFM3(2)
      write(6,*) 'SFM3(3) = ',SFM3(3)
      write(6,*) '==========================='
cccccccccc
C--
C--
C  Calculate sfermion masses.(to SSPTAB.inc)
C     SFM(I) = (SNL,SEL,SER,SUL,SUR,SDL,SDR)
C--
      IF(IS.NE.2) THEN
        CALL GTSFMS(AM0,AM2,AMU,TANB,ALF,S2W,ZM,SFM)
      ENDIF
        CALL GTSFMS3(AMU,TANB,ALF,S2W,ZM,FM3,ASFT,SFM,
     &            SFM3,SFM3M,SPHI,IS)
      CALL INOMIX3(AM,WM,S2W,ALF,AMU,M1,M2,TANB,SWM,SZM,
     .  OL,OR,ER,PHIL,PHIR,ON,UN,ETA,IFLG)
C      ENDIF
C--
C  CALCULATE HIGGS MASS
c--
      AT=ASFT(2)
      AB=ASFT(3)
      CALL MHIGGS(AMA, TANB, FM3, AT, AB, AMU,
     >            SFM3M, ALF, S2W, ZM, 
     >            ALPH, AMH1, AMH2,AMHC)
       SHM(1)=AMH1
       SHM(2)=AMH2
       SHM(3)=AMA
       SHM(4)=AMHC
       RETURN
       END  

 











