C***********************************************************************
C*
C*======================================----==
C* Subroutine PRUSER(IDREC,LEVEL,IDEBUG,IRET)
C*======================================----==
C*
C* (Function)
C*     User Event analysis routine.
C* (Input)
C*     IDREC   :  Record ID.
C*     ILEVEL  :  Level of analysis.
C*     IDEBUG  :  Debug flag.
C* (Output)
C*     IRET    :  Return flag.
C* (Update Record)
C*   95/02/12  K.Fujii  	Original version for debugging.
C*
C***********************************************************************
 
      SUBROUTINE PRUSER( IDREC, LEVEL, IDEBUG, IRET )
 
      CHARACTER*64    BKNAME(0:1)
      DATA BKNAME  / 'Spring:Parton_List',
     .		     'Generator:Particle_List' /

      PARAMETER      ( LNxDT = 100, MXxPT = 5000 )
      COMMON /SSWORK/ IELM(MXxPT,0:3), PV(4,500), RDATA(LNxDT,MXxPT)
      INTEGER*4       IDATA(LNxDT,MXxPT), NTK(0:3)
      EQUIVALENCE    (IDATA,RDATA)
      DATA IEVT /0/
C
C==================<<  Entry Point  >>================================
C
C--
C  Reset return code, Skip except EVENT record.
C--
      IRET = 0
      IF( IDREC .NE. 1 ) RETURN
      LEVEL  = LEVEL
      IDEBUG = IDEBUG
C--
      IEVT = IEVT + 1
      IF ( IEVT.GT.50 ) RETURN
C--
C  TBGET parton information.
C        RDATA(    1) = particle serial #
C             (    2) = particle ID
C             (    3) = mass
C             (    4) = charge
C             (    5) = Px
C             (    6) = Py
C             (    7) = Pz
C             (    8) = E
C             ( 9-11) = not used
C             (   12) = # daughter partons
C             (   13) = particle serial # of 1-st daughter
C             (   14) = particle serial # of parent
C             (15-16) = not used
C             (   17) = helicity
C             (   18) = color singlet group ID
C             (19-20) = not used
C--
      NPT = 0
      CALL UVZERO(2,NTK)
      DO 1000 IB = 1, 1
         CALL TBNOEL(1,BKNAME(IB),NE,IELM(1,IB))
         IF ( NE.LE.0 )                          GO TO 1000
         CALL UVZERO(4,PV(1,1))	 
         DO 100 IE = 1, NE
            NPT     = NPT + 1
            NTK(IB) = NTK(IB) + 1
            CALL TBGET(1,BKNAME(IB),IELM(IE,IB),NW,RDATA(1,NPT),IRT)
            CALL PLPRNT(IEVT,NPT,RDATA(1,NPT))
            IF ( INT(RDATA(12,NPT)).EQ.0 ) THEN
               CALL UADD4(RDATA(5,NPT),PV(1,1),PV(1,1))
	    ENDIF
100      CONTINUE
	 PRINT *, ' ----- '
	 PRINT *, ' PX, PY, PZ, E = ', PV(1,1), PV(2,1), PV(3,1), PV(4,1)
1000  CONTINUE
C--
C  That's it.
C--
      RETURN
      END
