      SUBROUTINE GETECN(NPT,ILP,LNxDT,RDATA,ECONTH,ECON,INCN)
 
      INTEGER*4  NPT, ILP, LNxDT, INCN(NPT)
      REAL   *4  RDATA(LNxDT,NPT), ECONTH, ECON
      DATA NCALL /0/
C
C==================<<  Entry Point  >>================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xPI    = ACOS(-1.)
         xRD2DG = 180/xPI
      ENDIF
C--
C  Sum up energies in the cone.
C--
      APLP   = UABSV3(RDATA(1,ILP))
      CSTHCT = COS(ECONTH/xRD2DG)
      CALL UVZERO(NPT,INCN)
      ECON   = 0
      DO 1000 IPT = 1, NPT
         IF ( IPT.EQ.ILP )                       GO TO 1000
         AP   = UABSV3(RDATA(1,IPT))
         CSTH = UDOT3(RDATA(1,ILP),RDATA(1,IPT))/APLP/AP
         IF ( CSTH.LT.CSTHCT )                   GO TO 1000
         INCN(IPT) = 1
         ECON      = ECON + RDATA(4,IPT)
1000  CONTINUE
C--
C  That's it.
C--
      RETURN
      END
