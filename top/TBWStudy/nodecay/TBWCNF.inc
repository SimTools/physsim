      PARAMETER       ( NEXLIN = 10 )
      COMMON /TBWCNF/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), DSGMDX, DPDEBM, BMENGY(0:4),
     .                 ISRBM
      REAL   *8        ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                 BMENGY
      INTEGER*4        IDPTCL, IHLCMB, ISRBM
