      PARAMETER       ( NEXLIN = 8 )
      COMMON /WWZCNF/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), DSGMDX, DPDEBM, BMENGY(0:4),
     .	               ISRBM, 
     .                 IMD1MN, IMD1MX, IMD2MN, IMD2MX, IMD3MN, IMD3MX,
     .                 ALAM, A, B, BTILDE
      REAL   *4        ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                 BMENGY,
     .                 ALAM, A, B, BTILDE
      INTEGER*4        IDPTCL, IHLCMB, ISRBM
     .                 IMD1MN, IMD1MX, IMD2MN, IMD2MX, IMD3MN, IMD3MX
