      PARAMETER       ( NEXLIN = 10 )
      COMMON /TTZCNF/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), DSGMDX, DPDEBM, BMENGY(0:4),
     .                 ISRBM, 
     .                 IMDZMN, IMDZMX, IMD1MN, IMD1MX, IMD2MN, IMD2MX,
     .                 POLP
      REAL   *4        ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                 BMENGY, POLP
      INTEGER*4        IDPTCL, IHLCMB, ISRBM,
     .                 IMDZMN, IMDZMX, IMD1MN, IMD1MX, IMD2MN, IMD2MX
