#ifndef WITH_DBD_STANDARD
      PARAMETER       ( NEXLIN = 10 )
#else
      PARAMETER       ( NEXLIN = 15)
#endif
      COMMON /TTBBCNF/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                 IDPTCL(3,NEXLIN), IHLCMB(NEXLIN),
     .                 PVOUT(0:3,NEXLIN), DSGMDX, DPDEBM, BMENGY(0:4),
     .                 ISRBM, 
     .                 IMD1MN, IMD1MX, IMD2MN, IMD2MX, POLP, USEALTCS
#ifdef WITH_DBD_STANDARD
     .               , Z3, Z4
      REAL*8           Z3, Z4
#endif
      REAL   *4        ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                 BMENGY, POLP
      INTEGER*4        IDPTCL, IHLCMB, ISRBM,
     .                 IMD1MN, IMD1MX, IMD2MN, IMD2MX,
     .                 USEALTCS
