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
      PARAMETER  ( MXxDAT = 81 )
      REAL   *8  BMDATA(2,0:MXxDAT)
C--
      DATA   SGEB / 0.005D0 /
      DATA ( BMDATA(I,0),I=1,2) /  2*0.D0 /
      DATA ((BMDATA(I,J),I=1,2),J= 1,20) /
     . 0.46917610D-03, 0.29999997D-01, 0.13529558D-02, 0.89999974D-01,
     . 0.26466998D-02, 0.14999998D+00, 0.42429529D-02, 0.20999998D+00,
     . 0.60631782D-02, 0.26999998D+00, 0.80721080D-02, 0.32999998D+00,
     . 0.10347933D-01, 0.38999999D+00, 0.12878112D-01, 0.44999999D+00,
     . 0.15549790D-01, 0.50999999D+00, 0.18368363D-01, 0.56999999D+00,
     . 0.21337438D-01, 0.63000000D+00, 0.24469260D-01, 0.69000000D+00,
     . 0.27896043D-01, 0.75000000D+00, 0.31481497D-01, 0.80999994D+00,
     . 0.35154924D-01, 0.86999995D+00, 0.39085470D-01, 0.92999995D+00,
     . 0.43139618D-01, 0.98999995D+00, 0.47354035D-01, 0.10499992D+01,
     . 0.51819582D-01, 0.11099997D+01, 0.56467917D-01, 0.11699991D+01/
      DATA ((BMDATA(I,J),I=1,2),J=21,40) /
     . 0.61408285D-01, 0.12299995D+01, 0.66563964D-01, 0.12900000D+01,
     . 0.72005093D-01, 0.13499994D+01, 0.77726722D-01, 0.14099998D+01,
     . 0.83691835D-01, 0.14699993D+01, 0.90102077D-01, 0.15299997D+01,
     . 0.96915781D-01, 0.15899992D+01, 0.10417134D+00, 0.16499996D+01,
     . 0.11192214D+00, 0.17099991D+01, 0.12023985D+00, 0.17699995D+01,
     . 0.12916076D+00, 0.18299999D+01, 0.13856155D+00, 0.18899994D+01,
     . 0.14857507D+00, 0.19499998D+01, 0.15900099D+00, 0.20099993D+01,
     . 0.17000854D+00, 0.20699997D+01, 0.18159306D+00, 0.21299992D+01,
     . 0.19364733D+00, 0.21899996D+01, 0.20635843D+00, 0.22500000D+01,
     . 0.21942562D+00, 0.23099995D+01, 0.23316675D+00, 0.23699999D+01/
      DATA ((BMDATA(I,J),I=1,2),J=41,60) /
     . 0.24709624D+00, 0.24299994D+01, 0.26161265D+00, 0.24899998D+01,
     . 0.27651101D+00, 0.25499992D+01, 0.29204160D+00, 0.26099997D+01,
     . 0.30791533D+00, 0.26699991D+01, 0.32411772D+00, 0.27299995D+01,
     . 0.34084654D+00, 0.27900000D+01, 0.35786092D+00, 0.28499994D+01,
     . 0.37542433D+00, 0.29099998D+01, 0.39330065D+00, 0.29699993D+01,
     . 0.40255284D+00, 0.30149994D+01, 0.41193247D+00, 0.30449991D+01,
     . 0.42150849D+00, 0.30749998D+01, 0.43127674D+00, 0.31049995D+01,
     . 0.44136357D+00, 0.31349993D+01, 0.45166171D+00, 0.31650000D+01,
     . 0.46204972D+00, 0.31949997D+01, 0.47280377D+00, 0.32249994D+01,
     . 0.48384172D+00, 0.32549992D+01, 0.49522501D+00, 0.32849998D+01/
      DATA ((BMDATA(I,J),I=1,2),J=61,80) /
     . 0.50713444D+00, 0.33149996D+01, 0.51948869D+00, 0.33449993D+01,
     . 0.53243959D+00, 0.33750000D+01, 0.54613018D+00, 0.34049997D+01,
     . 0.56078678D+00, 0.34349995D+01, 0.57676947D+00, 0.34649992D+01,
     . 0.59498650D+00, 0.34949999D+01, 0.62124234D+00, 0.35249996D+01,
     . 0.65583831D+00, 0.35549994D+01, 0.69752795D+00, 0.35849991D+01,
     . 0.74486202D+00, 0.36149998D+01, 0.79570866D+00, 0.36449995D+01,
     . 0.84795910D+00, 0.36749992D+01, 0.89336556D+00, 0.37049999D+01,
     . 0.92931241D+00, 0.37349997D+01, 0.95658273D+00, 0.37649994D+01,
     . 0.97634214D+00, 0.37949991D+01, 0.98944020D+00, 0.38249998D+01,
     . 0.99695528D+00, 0.38549995D+01, 0.99987316D+00, 0.38849993D+01/
      DATA ((BMDATA(I,J),I=1,2),J=81,81) /
     . 0.99999988D+00, 0.39150000D+01/
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
