      PARAMETER       ( NEXLIN = 10 )
      COMMON /XNXNCN/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), DSGMDX, DPDEBM, BMENGY(0:4),
     .		       ISRBM, IMD1MN, IMD1MX, IMD2MN, IMD2MX
      REAL   *8        ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                 BMENGY
      INTEGER*4        IDPTCL, IHLCMB, ISRBM,
     .                 IMD1MN, IMD1MX, IMD2MN, IMD2MX
