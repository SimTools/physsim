      PARAMETER       ( NEXLIN = 10 )
      COMMON /SFSFCN/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), DSGMDX, DPDEBM, BMENGY(0:4),
     .                 HTAUM, HTAUP,
     .                 IHANDM, IHANDP, IGENSF, IDOTAU, ISRBM
      REAL   *4        ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                 BMENGY, HTAUM, HTAUP
      INTEGER*4        IDPTCL, IHLCMB, IHANDM, IHANDP, IGENSF, 
     .                 IDOTAU, ISRBM
