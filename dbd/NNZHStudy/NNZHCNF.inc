      PARAMETER       ( NEXLIN = 7 )
C .. NEXLIN=8, 9 for e+/e- after beamstrahlung
C .. NEXLIN=10, 11 for ISR gamma
      COMMON /NNZHCNF/  ROOTS, POLE, XPHASE(3,NEXLIN-3),
     .                  PVOUT(0:3,NEXLIN+4), DSGMDX, DPDEBM,
     .                  BMENGY(0:4), POLP,
#ifdef WITH_DBD_STANDARD
     .                  Z3, Z4,
#endif
     .                  IDPTCL(3,NEXLIN), IHLCMB(NEXLIN), ISRBM,
     .                  IMDMN, IMDMX
#ifdef WITH_DBD_STANDARD
      REAL*8            Z3, Z4
#endif
      REAL   *8         ROOTS, POLE, XPHASE, PVOUT, DSGMDX, DPDEBM,
     .                  BMENGY, POLP
      INTEGER*4         IDPTCL, IHLCMB, ISRBM,
     .                  IMDMN, IMDMX
