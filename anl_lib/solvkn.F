      SUBROUTINE SOLVKN(RS,AMP,AMD,P1,P2,EN,IRET)
 
      REAL   *4  RS, AMP, AMD, P1(4), P2(4), EN(3,2)
      INTEGER*4  IRET
C--
      REAL   *4  EN1(3), EN2(3), EE(3,3), ENEE(3,2)
C
C==================<<  Entry Point  >>================================
C
C--
C  Initialization.
C--
      IRET = -1
C--
C  Decide the parent direction of P1.
C--
      CALL UNRMV3(P1,EN1)
      CALL UNRMV3(P2,EN2)
      BT     =  SQRT( (1-2*AMP/RS)*(1+2*AMP/RS) )
      E1     =  P1(4)
      E2     =  P2(4)
      ENEN1  =  ( RS*E1 - (AMP-AMD)*(AMP+AMD) )/E1/RS/BT
      ENEN2  = -( RS*E2 - (AMP-AMD)*(AMP+AMD) )/E2/RS/BT
      EN1EN2 =  UDOT3(EN1,EN2)
C--
      AEE2    =  (1-EN1EN2)*(1+EN1EN2)
      IF ( AEE2.LE.0. )                          GO TO 9999
      AEE2    =  1/SQRT(AEE2)
C--
      EE(1,1) = EN1(1)
      EE(2,1) = EN1(2)
      EE(3,1) = EN1(3)
      EE(1,2) = AEE2*( EN2(1) - EN1EN2*EN1(1) )
      EE(2,2) = AEE2*( EN2(2) - EN1EN2*EN1(2) )
      EE(3,2) = AEE2*( EN2(3) - EN1EN2*EN1(3) )
      CALL UNCROS(EE(1,1),EE(1,2),EE(1,3))
C--
      ENEE(1,1)  =  ENEN1
      ENEE(2,1)  =  AEE2*( ENEN2 - EN1EN2*ENEN1 )
      ENEE(3,1)  =  SQRT(MAX( 1 - ENEE(1,1)**2 - ENEE(2,1)**2, 0.) )
      ENEE(1,2)  =  ENEE(1,1)
      ENEE(2,2)  =  ENEE(2,1)
      ENEE(3,2)  = -ENEE(3,1)
C--
      DO 100 K = 1, 2
         DO 10 I = 1, 3
            EN(I,K) = 0
            DO 1 J = 1, 3
               EN(I,K) = EN(I,K) + ENEE(J,K)*EE(I,J)
1           CONTINUE
10       CONTINUE
100   CONTINUE
C--
C  That's it.
C--
      IRET = 0
9999  RETURN
      END
