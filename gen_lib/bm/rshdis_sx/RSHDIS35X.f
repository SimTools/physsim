C
C  *INCLUDE UDSRCH
C        WRITE(21,'(''NEWFRAME'')')
C        WRITE(21,'(''SET ORDER X Y'')')
C        DO 10 I = 1, 100
C           X = 0.01*I
C           CALL RSHDIS(X,1,Y)
C           WRITE(21,'(2E15.8)') Y, X
C  10    CONTINUE
C        WRITE(21,'(''JOIN'')')
C        END
C*
C* This version for JLC-I X-band(150) parameter as of 93/04/18
C* with the bremsstrahlung:
C*                   _
C*                  / 1
C*    sig_RC(s) =  /   dz b [ z**(b-1)(1+3b/4) - (1-z/2) ] sig((1-z)s)
C*               _/ 0
C* The beam energy spread is assumed to be flat and has FWHM of 1.0%.
C*
      SUBROUTINE RSHDIS(X,ITYP,Y)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      INTEGER*4  ITYP
      REAL   *4  X, Y
C--
      DATA   SGEB / 0.0035D0 /
      PARAMETER  ( MXxDAT = 94 )
      REAL   *8  BMDATA(2,0:MXxDAT)
      DATA ( BMDATA(I,0),I=1,2) /  2*0.D0 /
      DATA ((BMDATA(I,J),I=1,2),J=  1, 20) /
     .  .34437698D-03,  .29999999D-01,  .11870501D-02,  .90000004D-01,
     .  .24122407D-02,  .15000001D+00,  .38882925D-02,  .20999999D+00,
     .  .56198803D-02,  .27000001D+00,  .77339043D-02,  .33000001D+00,
     .  .10082947D-01,  .38999999D+00,  .12487216D-01,  .44999999D+00,
     .  .15171436D-01,  .50999999D+00,  .17737417D-01,  .56999999D+00,
     .  .20741871D-01,  .63000000D+00,  .23772009D-01,  .69000000D+00,
     .  .27063387D-01,  .75000000D+00,  .30548772D-01,  .81000000D+00,
     .  .33976872D-01,  .87000000D+00,  .37889518D-01,  .93000001D+00,
     .  .41880738D-01,  .99000001D+00,  .45855127D-01,  .10500000D+01,
     .  .50094124D-01,  .11100000D+01,  .54635599D-01,  .11700000D+01/
      DATA ((BMDATA(I,J),I=1,2),J= 21, 40) /
     .  .59222151D-01,  .12300000D+01,  .64188905D-01,  .12900000D+01,
     .  .69216684D-01,  .13500000D+01,  .74455276D-01,  .14100000D+01,
     .  .80121085D-01,  .14700000D+01,  .86064570D-01,  .15300000D+01,
     .  .92385612D-01,  .15900000D+01,  .99073082D-01,  .16500000D+01,
     .  .10629798D+00,  .17100000D+01,  .11389209D+00,  .17700000D+01,
     .  .12194344D+00,  .18300000D+01,  .13052420D+00,  .18900000D+01,
     .  .13958256D+00,  .19500000D+01,  .14878361D+00,  .20100000D+01,
     .  .15876825D+00,  .20699999D+01,  .16907871D+00,  .21300001D+01,
     .  .17991453D+00,  .21900001D+01,  .19103776D+00,  .22500000D+01,
     .  .20232174D+00,  .23099999D+01,  .21426715D+00,  .23699999D+01/
      DATA ((BMDATA(I,J),I=1,2),J= 41, 60) /
     .  .22675936D+00,  .24300001D+01,  .23933053D+00,  .24900000D+01,
     .  .25235796D+00,  .25500000D+01,  .26536310D+00,  .26099999D+01,
     .  .27914733D+00,  .26700001D+01,  .29283342D+00,  .27300000D+01,
     .  .30698615D+00,  .27900000D+01,  .32134384D+00,  .28499999D+01,
     .  .33568063D+00,  .29100001D+01,  .34999892D+00,  .29700000D+01,
     .  .35764167D+00,  .30150001D+01,  .36529404D+00,  .30450001D+01,
     .  .37321264D+00,  .30750000D+01,  .38063285D+00,  .31050000D+01,
     .  .38837019D+00,  .31350000D+01,  .39620182D+00,  .31650000D+01,
     .  .40404463D+00,  .31949999D+01,  .41207501D+00,  .32249999D+01,
     .  .42028633D+00,  .32550001D+01,  .42820507D+00,  .32850001D+01/
      DATA ((BMDATA(I,J),I=1,2),J= 61, 80) /
     .  .43632859D+00,  .33150001D+01,  .44462192D+00,  .33450000D+01,
     .  .45320043D+00,  .33750000D+01,  .46179986D+00,  .34050000D+01,
     .  .47028947D+00,  .34349999D+01,  .47917986D+00,  .34649999D+01,
     .  .48828349D+00,  .34949999D+01,  .49738002D+00,  .35250001D+01,
     .  .50687194D+00,  .35550001D+01,  .51650238D+00,  .35850000D+01,
     .  .52642602D+00,  .36150000D+01,  .53667957D+00,  .36450000D+01,
     .  .54733056D+00,  .36750000D+01,  .55827534D+00,  .37049999D+01,
     .  .56952560D+00,  .37349999D+01,  .58155823D+00,  .37650001D+01,
     .  .59477973D+00,  .37950001D+01,  .60889781D+00,  .38250000D+01,
     .  .62595272D+00,  .38550000D+01,  .65092939D+00,  .38850000D+01/
      DATA ((BMDATA(I,J),I=1,2),J= 81, 94) /
     .  .68369466D+00,  .39150000D+01,  .72326410D+00,  .39449999D+01,
     .  .76755601D+00,  .39749999D+01,  .81535685D+00,  .40050001D+01,
     .  .86374807D+00,  .40349998D+01,  .90482259D+00,  .40650001D+01,
     .  .93739253D+00,  .40949998D+01,  .96222526D+00,  .41250000D+01,
     .  .97985315D+00,  .41550002D+01,  .99154538D+00,  .41849999D+01,
     .  .99782503D+00,  .42150002D+01,  .99997646D+00,  .42449999D+01,
     .  .10000000D+01,  .42750001D+01,  .10000000D+01,  .43049998D+01/
C
C========< Entry Point >================================================
C
C--
C  Decide which bin X belongs to.
C--
      CALL UDSRCH(MXxDAT,2,1,BMDATA(1,1),DBLE(X),IP)
C--
C  Calculate Y.
C     FRAC = MAX(1.E-6,1+5*SGEB-FRAC)
C     FRAC = -LOG(FRAC)
C--
      IF ( IP.NE.0 .AND. IP.NE.MXxDAT ) THEN
         XL  = BMDATA(1,IP)
         DX  = BMDATA(1,IP+1) - XL
         DF  = (BMDATA(2,IP+1)-BMDATA(2,IP))/DX
         VAL = BMDATA(2,IP) + (X-XL)*DF
      ELSE
         VAL = BMDATA(2,IP)
      ENDIF
      Y   = 1 + 5*SGEB - EXP(-VAL)
C--
C  That's it.
C--
      RETURN
      END
