CC********************************************************************CC
C*                                                                    *C
C*=================================-----===                           *C
C* Subroutine LAMBDA(ALFS,ENF,TYPE,QOVLM)                             *C
C*=================================-----===                           *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    Calculates Q/LAM for given alpha_s with various formulae.       *C
C* (Inputs)                                                           *C
C*    ALFS    :(R*4);   alpha-s.                                      *C
C*    ENF     :(R*4);   # active quark flavors.                       *C
C*    TYPE    :(C*4);   formula type.                                 *C
C*                  ;    'RM03' = 3-loop by Marshall.                 *C
C*                  ;    'KF03' = 3-loop by Marshall with usual       *C
C*                  ;             boundary condition.                 *C
C*                  ;    'SG03' = 3-loop by Gorishny.                 *C
C*                  ;    'TT03' = 3-loop by Tauchi.                   *C
C*                  ;    'TT02' = 2-loop by PDG.                      *C
C*                  ;    'MK02' = 2-loop by MK-II.                    *C
C* (Output)                                                           *C
C*    QOVLM   :(R*4);   Q/Lambda_MS_bar.                              *C
C* (Relation)                                                         *C
C*    Invokes no subroutine or function.                              *C
C* (Update Record)                                                    *C
C*    8/09/89  K.Fujii          Original version                      *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE LAMBDA(ALFS,ENF,TYPE,QOVLM)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*4       ALFS, ENF, QOVLM
      CHARACTER*4  TYPE
      DATA NCALL   / 0 /
      DATA TEST    / 1.D-5 /
      DATA MXxTRY  / 500 /
C
C========< Entry Point >================================================
C
C--
C  Prepare some numerical constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         xPI   = ACOS(-1.D0)
         xHPI  = xPI/2
         x2PI  = 2*xPI
         x4PI  = 4*xPI
      ENDIF
C--
C  Prepare coefficents for beta-function.
C     BT0 = beta-0  ( 1-loop; sheme independent   )
C     BT1 = beta-1  ( 2-loop; sheme independent   )
C     BT2 = beta-2  ( 3-loop; Lambda-MS-bar sheme )
C--
      BT0 = 11 - 2*ENF/3
      BT1 = 2*( 51 - 19*ENF/3 )
      BT2 = ( 2857 - 5033*ENF/9 + 325*ENF**2/27 )/2
C--
C  Branch on types of alpha-s formulae.
C--
      X   = x4PI/(BT0*ALFS)
C--
C  Marshall's formula.
C--
      IF ( TYPE.EQ.'RM03' .OR. TYPE.EQ.'KF03' .OR. TYPE.EQ.'SG03' ) THEN
         RDLT = 1/SQRT( (4*BT2*BT0-BT1**2)/BT0**4 )
         Y    = X - 0.5*(BT1/BT0**2)
     .              *LOG( X**2 + (BT1/BT0**2)*X + BT2/BT0**3 )
     .            - RDLT*( (BT1**2-2*BT2*BT0)/BT0**4 )
     .              *ATAN( RDLT*( BT1/BT0**2 + 2*(BT2/BT0**3)/X ) )
         IF ( TYPE.EQ.'KF03' ) THEN
            Y = Y + 0.5*(BT1/BT0**2)*LOG(BT2/BT0**3)
     .            + xHPI*RDLT*( (BT1**2-2*BT2*BT0)/BT0**4 )
         ELSE IF ( TYPE.EQ.'SG03' ) THEN
            Y = Y + RDLT*( (BT1**2-2*BT2*BT0)/BT0**4 )
     .              *ATAN( RDLT*BT1**2/BT0**4 )
     .            - (BT1/BT0**2)*LOG(BT1/BT0**2)
         ENDIF
         QOVLM = EXP(Y/2)
C--
C  Tauchi's formula.
C--
      ELSE IF ( TYPE.EQ.'TT03' .OR. TYPE.EQ.'TT02' ) THEN
         B1B02 = BT1/BT0**2
         B12   = BT1**2
         B0B2  = BT0*BT2
         B04   = 1/BT0**4
C--
         Z     = 1.D-1
         RX    = 1/X
         NTRY  = 0
         ADFS  = 1.D20
C--
1        NTRY  = NTRY + 1
         ALZ   = LOG(Z)
         F     = Z + B1B02*Z**2*ALZ
         DFDZ  = 1 + B1B02*(2*Z*ALZ+1/Z)
         IF ( TYPE.EQ.'TT03' ) THEN
            F     = F + (B12*(ALZ**2+ALZ-1)+B0B2)*B04*Z**3
            DFDZ  = DFDZ + B12*(2*ALZ/Z+1/Z)*B04*Z**3
     .                   + (B12*(ALZ**2+ALZ-1)+B0B2)*3*B04*Z**2
         ENDIF
         DF = RX - F
         IF ( NTRY.GT.MXxTRY ) THEN
            PRINT *, ' NTRY over '
            STOP
         ENDIF
         IF ( ABS(DF).LT.TEST )                  GO TO 10
         DZ   = DF/DFDZ
         Z    = Z + DZ
         ADFS = ABS(DF)
                                                 GO TO 1
10       QOVLM = EXP(1/(2*Z))
C        PRINT *, ' NTRY = ', NTRY
C--
C  MARK-II formula.
C--
      ELSE IF ( TYPE.EQ.'MK02' ) THEN
         B1B02 = BT1/BT0**2
C--
         Y     = 1.D-1
         NTRY  = 0
         ADFS  = 1.D20
C--
2        NTRY  = NTRY + 1
         F     = Y + B1B02*LOG(Y)
         DFDY  = 1 + B1B02/Y
         DF    = X - F
         IF ( NTRY.GT.MXxTRY ) THEN
            PRINT *, ' NTRY over '
            STOP
         ENDIF
         IF ( ABS(DF).LT.TEST )                  GO TO 20
         DY   = DF/DFDY
         Y    = Y + DY
         ADFS = ABS(DF)
                                                 GO TO 2
20       QOVLM = EXP(Y/2)
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
