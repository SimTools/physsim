      SUBROUTINE IDK2PD(IDKF,IDPD)
 
      INTEGER*4  IDKF(3), IDPD
      INTEGER*4  IDPDG(3,2,2)
      DATA IDPDG / 12, 14, 16,    11, 13, 15,
     .              2,  4,  6,     1,  3,  5 /
 
      IG = IDKF(1)
      IT = IDKF(2)
      LQ = IDKF(3)
      IDPD = IDPDG(IG,IT,LQ)
      RETURN
      END
