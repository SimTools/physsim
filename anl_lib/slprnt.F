      SUBROUTINE SLPRNT(IEVT,IPT,RDATA)
 
      INTEGER*4  IEVT, IPT
      REAL   *4  RDATA(100)
      DATA IEVTS / -999 /
C
C==================<<  Entry Point  >>================================
C
C--
C  Header.
C--
      IF ( IEVT.NE.IEVTS ) THEN
         IEVTS = IEVT
         WRITE(*,'(1H0,''********* IEVT = '',I8,'' *************'',
     .                 '' Spring:Parton_List *******************'')')
     .                    IEVT
         WRITE(*,'(/,'' IPT    P_ID      Mass      Q       Px       Py''
     .               ''       Pz        E       Vx       Vy       Vz'',
     .                 '' NDT ID1 IPR  Helic.  IC  IS'')')
         WRITE(*,'('' '',120(''-''))')
      ENDIF
C--
C  Print contents of 'Spring:Particle_List'.
C--
      IPID = RDATA( 2)
      Q    = RDATA( 4)
      NDT  = RDATA(12)
      ID1  = RDATA(13)
      IPR  = RDATA(14)
      ICG  = RDATA(18)
      ISG  = RDATA(19)
      WRITE(*,'(1I4,1I8,1F10.4,1F7.3,7F9.3,3I4,F8.4,2I4)')
     .    IPT, IPID, RDATA(3), Q, (RDATA(K),K=5,11), NDT, ID1, IPR,
     .    RDATA(17), ICG, ISG
C--
C  That's it.
C--
      RETURN
      END
