CC********************************************************************CC
C*                                                                    *C
C*==================================----===                           *C
C* Subroutine ALPHAS(QOVLM,ENF,TYPE,ALFS)                             *C
C*==================================----===                           *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    Calculates alpha_s for given Q/LAM with various formulae.       *C
C* (Inputs)                                                           *C
C*    QOVLM   :(R*4);   Q/Lambda_MS_bar.                              *C
C*    ENF     :(R*4);   # active quark flavors.                       *C
C*    TYPE    :(C*4);   formula type.                                 *C
C*                  ;    'TT03' = 3-loop by Tauchi.                   *C
C*                  ;    'TT02' = 2-loop by PDG.                      *C
C*                  ;    'TT01' = 1-loop by Tauchi.                   *C
C*                  ;    'MK02' = 2-loop by MK-II.                    *C
C* (Output)                                                           *C
C*    ALFS    :(R*4);   alpha-s.                                      *C
C* (Relation)                                                         *C
C*    Invokes no subroutine or function.                              *C
C* (Update Record)                                                    *C
C*    8/09/89  K.Fujii          Original version                      *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE ALPHAS(QOVLM,ENF,TYPE,ALFS)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*4       QOVLM, ENF, ALFS
      CHARACTER*4  TYPE
      DATA NCALL   / 0 /
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
C     BT0 = beta-0  ( 1-loop; scheme independent  )
C     BT1 = beta-1  ( 2-loop; scheme independent  )
C     BT2 = beta-2  ( 3-loop; Lambda-MS-bar sheme )
C--
      BT0 = 11 - 2*ENF/3
      BT1 = 2*( 51 - 19*ENF/3 )
      BT2 = ( 2857 - 5033*ENF/9 + 325*ENF**2/27 )/2
C--
C  Branch on types of alpha-s formulae.
C--
      IF ( TYPE.EQ.'MK02' ) THEN
         X = 2*LOG(QOVLM) + (BT1/BT0**2)*LOG(2*LOG(QOVLM))
         X = 1/(BT0*X)
      ELSE
         Y = 1/(BT0*2*LOG(QOVLM))
         Z = LOG(2*LOG(QOVLM))
         X = Y
         IF ( TYPE.EQ.'TT02' .OR. TYPE.EQ.'TT03' ) THEN
            X = X - (BT1*Z/BT0)*Y**2
            IF ( TYPE.EQ.'TT03' ) THEN
               X = X + ( (Z**2-Z-1)*(BT1/BT0)**2 + BT2/BT0 )*Y**3
            ENDIF
         ENDIF
      ENDIF
      ALFS = x4PI*X
C--
C  That's it.
C--
      RETURN
      END
