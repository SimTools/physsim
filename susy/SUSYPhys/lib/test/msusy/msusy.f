CINCLUDE SFMASS
CINCLUDE XCMASS
CINCLUDE XNMASS
CINCLUDE USORTD
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
 
      REAL*8      AM0, AM2, AMU, BT, AMW, AMZ, ALF, S2W, AMSF2(7),
     .            AMXC(2), PHIL, PHIR, EPSR,
     .            AMXN(4), ON(4,4)
      COMPLEX*16  ETA(4)
 
      AM0 = 170
      AM2 = 200
      AMU = 400
      TNB = +2
      BT  = ATAN(TNB)
      ALF = 1/128.D0
      AMZ = 91.17D0
      S2W = 0.230D0
      AMW = AMZ*SQRT(1-S2W)
      CALL SFMASS(AM0,AM2,AMU,BT,ALF,S2W,AMZ,AMSF2)
      CALL XCMASS(AM2,AMU,BT,AMW,AMXC,PHIL,PHIR,EPSR)
      CALL XNMASS(AM2,AMU,BT,S2W,AMZ,AMXN,ON,ETA)
      PRINT *, 'AM0 = ', AM0
      PRINT *, 'AM2 = ', AM2
      PRINT *, 'AMU = ', AMU
      PRINT *, 'TNB = ', TNB
      PRINT *, 'SWM(1) = ', AMXC(1)
      PRINT *, '   (2) = ', AMXC(2)
      PRINT *, 'SZM(1) = ', AMXN(1)
      PRINT *, '   (2) = ', AMXN(2)
      PRINT *, '   (3) = ', AMXN(3)
      PRINT *, '   (4) = ', AMXN(4)
      PRINT *, 'SLR    = ', SQRT(AMSF2(7))
      PRINT *, 'SLL    = ', SQRT(AMSF2(5))
      PRINT *, 'SNL    = ', SQRT(AMSF2(6))
      END
