      PARAMETER       ( NEXLIN = 6 )
      COMMON /NNZZCNF/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), 
     .                 ISRBM, IMD1MN, IMD1MX, IMD2MN, IMD2MX
      REAL   *8        ROOTS, POLE, XPHASE, PVOUT
      INTEGER*4        IDPTCL, IHLCMB, ISRBM, 
     .                 IMD1MN, IMD1MX, IMD2MN, IMD2MX
