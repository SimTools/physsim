C----------------------------------------------------------------------C
      SUBROUTINE HGSMIX( SIN2W, WM, TANB, PSMAS, HG0MAS, HGCMAS, TANA)
C----------------------------------------------------------------------C
C  purpose : give information about neutral/charged Higgs bosons
C    ( inputs & outputs are real*4 )
C  inputs  : SIN2W / weak mixing angle
C            WM    / W-boson mass
C            PSMAS / pseudoscalar mass
C            TANB  / ratio of vev. ( v2/v1 )
C  outputs : HG0MAS(i) i=1..2 / neutral Higgs mass (i=1:lighter )
C            HGCMAS           / charged Higgs mass( physical )
C            TANA   / tan of Higgs boson mixing angle
C----------------------------------------------------------------------C
      REAL*4 HG0MAS(2), HG0MS2(2)
C     DATA IIII/0/
C----------------------------------------------------------------------C
      COS2W = 1.-SIN2W
      COSW  = SQRT(COS2W)
      ZM    = WM/COSW
      ZM2   = ZM**2
      PSMAS2= PSMAS**2
 
C     BETA  = DATAN2( DABS(DBLE(TANB)), DSIGN(1.D0,DBLE(TANB)) )
C     COSB  = COS(BETA)
C     SINB  = SIN(BETA)
C     COSDB = (COSB-SINB)*(COSB+SINB)
C     SINDB = 2.*SINB*COSB
C     SIN2B = SINB**2
C     COS2B = COSB**2
      TAN2B = TANB**2
      COS2B = 1./(1.+TAN2B)
      COSB  = SIGN( SQRT(COS2B), TANB )
      SIN2B = 1.-COS2B
      SINB  = SQRT(SIN2B)
      COSDB = (COSB-SINB)*(COSB+SINB)
      SINDB = 2.*SINB*COSB
 
      HGCMAS = SQRT( PSMAS2 + WM**2 )
      AMS    = PSMAS2 + ZM2
      DELTM  = SQRT( ((PSMAS+ZM)*(PSMAS-ZM))**2
     &                + 4.*PSMAS2 * ZM2 *SINDB**2 )
 
      HG0MS2(1) =  0.5*( AMS - DELTM )
      HG0MS2(2) =  0.5*( AMS + DELTM )
      DO 10 I=1,2
10    HG0MAS(I) = SQRT( HG0MS2(I) )
 
CCCC
C     IF( IIII.EQ.0 ) THEN
C     WRITE(6,*) 'COSB=',COSB
C     WRITE(6,*) 'SINB=',SINB
C     WRITE(6,*) 'ZM=',ZM
C     WRITE(6,*) 'PSMAS=',PSMAS
C     WRITE(6,*) 'DELTM=',DELTM
C     WRITE(6,*) 'H0MAS=',HG0MAS
C     WRITE(6,*) 'HCMAS=',HGCMAS
C     IIII=1
C     ENDIF
CCCC
 
      TANA = ( HG0MS2(1) -PSMAS2*COS2B -ZM2*SIN2B )
     &      / ( PSMAS2+ZM2 )/SINB/COSB
C     ALPHA = ATAN2( ABS(TANA), SIGN(1.,TANA) )
      ALPHA = DATAN2( DABS(DBLE(TANA)), DSIGN(1.D0,DBLE(TANA)) )
 
 
      RETURN
      END
