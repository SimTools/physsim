C *INCLUDE  UDSRCH
C      IMPLICIT REAL*8 (A-H,O-Z)
C      DO 100 I = 1, 100
C        X = 0.01*FLOAT(I)
CCC    X = 0.5
C      CALL EPROBX(X,11, XBEAM)
C
C       PRINT *,' X=',X,' XBEAM =',XBEAM
C100   CONTINUE
C      STOP
C      END
C**********************************************************************
C*
C*  -----------------------------=====
C*   SUBROUTINE  EPROBX(X, NTYP, XBEAM)
C*  -----------------------------=====
C*(Function)
C*   Generate Beam energy
C*   E at the collision = Ebeam*Xbeam
C*
C*(Input)
C*   X  : Integral variable, 0 < X < 1
C*   NTYP  : 1 = Ebeam = 250GeV
C*           2 =       = 500GeV
C*           3 =       = 750GeV
C*          11 = Ebeam = 250GeV with 0.5% initial beam energy spread.
C*          12 =       = 500GeV           (Gaussian)
C*          13 =       = 750GeV
C*          21 = Ebeam = 250GeV with 1%   initial beam energy spread.
C*          22 =       = 500GeV           (Full width, Sqauare)
C*          23 =       = 750GeV
C*          31 = JLC-I Sband 150 GeV, 0.01% width
C*          32 =       Cband   
C*          33 =       Xband 
C*          34 = JLC-I Sband 250 GeV, 0.01% width
C*          35 =       Cband   
C*          36 =       Xband 
C*          41 = JLC-I Sband 150 GeV, 0.5% width
C*          42 =       Cband   
C*          43 =       Xband 
C*          44 = JLC-I Sband 250 GeV, 0.5% width
C*          45 =       Cband   
C*          46 =       Xband 
C*(Output)
C*   EPROBX  : function to generate beam spectrum.
C*             0 < EPROBX < 1.025 ( Corresponds to 5 sigma when
C*                                  sigma is 0.5%.)
C*(Required subroutine)
C*   Sub. UDSRCH ; Do binary search of real*8 variable.
C*
C*(Author)
C*   A. Miyamoto  16-Oct-1990  Original version.
C*                11-Jun-1991  For JLC parameter on May 1991.
C*                31-Jul-1991  Add an option, 1% square Ebeam spread
C*                 6-Feb-1993  Add JLC-I parameters
C*
C**********************************************************************
C*
      SUBROUTINE EPROBX( X, NTYP, XBEAM )
C*
      IMPLICIT REAL*8 (A-H,O-Z)
C****************************************************
C*  Block data to generate Beamstraulung spectrum
C*****************************************************
      PARAMETER (NO$DAT=0200)
      REAL*8     XB(2,NO$DAT,20)
C
C ... Type = 1
      INCLUDE 'beaminc/eprob.blk.sd250_1D4'
C ... Type = 2
      INCLUDE 'beaminc/eprob.blk.sd500_1D4'
C ... Type = 3
      INCLUDE 'beaminc/eprob.blk.sd750_1D4'
C
C*****************************************************
C*  Block data to generate Beamstraulung spectrum
C*  250 GEV, 1% WIDTH
C*****************************************************
C ... Type = 4
      INCLUDE 'beaminc/eprob.blk.type4'
C ... Type = 5
      INCLUDE 'beaminc/eprob.blk.type5'
C ... Type = 6
      INCLUDE 'beaminc/eprob.blk.type6'
C
C*****************************************************
C   Beam Parameter for JLC-I   0.01% width
C****************************************************
C ... Type = 7
      INCLUDE 'beaminc/eprob.blk.s150_1D4'
C ... Type = 8
      INCLUDE 'beaminc/eprob.blk.c150_1D4'
C ... Type = 9
      INCLUDE 'beaminc/eprob.blk.x150_1D4'

C ... Type = 10
      INCLUDE 'beaminc/eprob.blk.s250_1D4'
C ... Type = 11
      INCLUDE 'beaminc/eprob.blk.c250_1D4'
C ... Type = 12
      INCLUDE 'beaminc/eprob.blk.x250_1D4'

C*****************************************************
C   Beam Parameter for JLC-I   0.5% width
C****************************************************
C ... Type = 13
      INCLUDE 'beaminc/eprob.blk.s150_5D3'
C ... Type = 14
      INCLUDE 'beaminc/eprob.blk.c150_5D3'
C ... Type = 15
      INCLUDE 'beaminc/eprob.blk.x150_5D3'

C ... Type = 16
      INCLUDE 'beaminc/eprob.blk.s250_5D3'
C ... Type = 17
      INCLUDE 'beaminc/eprob.blk.c250_5D3'
C ... Type = 18
      INCLUDE 'beaminc/eprob.blk.x250_5D3'
 
C *********************************************
C  For Spectrum with 0.5% gaussian spread.
C *********************************************
      PARAMETER (IDBMAX= 3)
      REAL*8 XX(100,IDBMAX), FX(100,IDBMAX)
 
      INCLUDE 'beaminc/sd250b.f'
      INCLUDE 'beaminc/sd500b.f'
      INCLUDE 'beaminc/sd750b.f'
C
C*
C==========< Entry Point >=============================================
C*
C
      IF( NTYP .GT. 10 .AND. NTYP .LT. 20   ) GO TO 200
C*
C*
      IF(NTYP.LT.11) THEN
         ITYP = NTYP
         XFACT=1.D0
      ELSEIF(NTYP.LT.31) THEN
         ITYP = NTYP - 20 + 3
         XFACT=1.005D0
      ELSEIF(NTYP.LT.41) THEN
         XFACT=1.D0
         ITYP = NTYP - 30 + 6
      ELSEIF(NTYP.LT.51) THEN
         ITYP = NTYP - 40 + 12
         XFACT=1.005D0
      ELSE
         PRINT *,'Parameter error in sub. EPROBX  beam type was =',NTYP
      ENDIF

      CALL UDSRCH( NO$DAT, 2, 2, XB(1,1,ITYP), X, NPNT )
C
      IF( NPNT .LE. 0 ) THEN
        XBEAM = 0.0
      ELSEIF( NPNT .GE. NO$DAT) THEN
        XBEAM = 1.0
      ELSE
        XLOG =  XB(1,NPNT,ITYP)
     >       + (XB(1,NPNT+1,ITYP)-XB(1,NPNT,ITYP))
     >       / (XB(2,NPNT+1,ITYP)-XB(2,NPNT,ITYP))
     >       * (X - XB(2,NPNT,ITYP))
CC      PRINT *,' XLOG '
        XBEAM = 1.D0 - 10.D0**(-1.*XLOG)
        IF( XBEAM .GT. 1.0D0 ) XBEAM = 1.0D0
        IF( XBEAM .LT. 0.0D0 ) XBEAM = 0.0D0
      ENDIF
C ... Adjust beam energy
      XBEAM = XBEAM*XFACT
 
      RETURN
C
C For spectrum with 0.5% gaussian spread.
C
200   CONTINUE
      INDATA = 100
      ITYPE = MOD(NTYP,10)
      CALL UDSRCH( INDATA, 1, 1, FX(1, ITYPE), X, NPNT )
CC    PRINT *,' NPNT =',NPNT
C
      IF( NPNT .LE. 0 ) THEN
        XBEAM = 0.0
      ELSEIF( NPNT .GE. INDATA) THEN
        XBEAM = XX(100,ITYPE)
      ELSE
        XBEAM =  XX(NPNT,ITYPE)
     >       + (XX(NPNT+1,ITYPE)-XX(NPNT,ITYPE))
     >       / (FX(NPNT+1,ITYPE)-FX(NPNT,ITYPE))
     >       * (X - FX(NPNT,ITYPE))
        IF( XBEAM .GT. XX(100,ITYPE) ) XBEAM = XX(100,ITYPE)
        IF( XBEAM .LT. 0.0D0 ) XBEAM = 0.0D0
      ENDIF
 
      RETURN
      END
