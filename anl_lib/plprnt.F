      SUBROUTINE PLPRNT(IEVT,IPT,RDATA)
 
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
     .                 '' Generator:Particle_List **************'')')
     .                    IEVT
         WRITE(*,'(/,'' IPT    P_ID      MASS   Q       PX       PY''
     .               ''       Pz        E       Vx       Vy       Vz'',
     .                 '' NDT ID1 IPR      CT      DL'')')
         WRITE(*,'('' '',120(''-''))')
      ENDIF
C--
C  Print contents of 'Generator:Particle_List'.
C--
      IPID = RDATA( 2)
      IQ   = RDATA( 4)
      NDT  = RDATA(12)
      ID1  = RDATA(13)
      IPR  = RDATA(14)
      WRITE(*,'(1I4,1I8,1F10.4,1I4,7F9.3,3I4,2F8.4)')
     .    IPT, IPID, RDATA(3), IQ, (RDATA(K),K=5,11), NDT, ID1, IPR,
     .    (RDATA(K),K=16,17)
C--
C  That's it.
C--
      RETURN
      END
