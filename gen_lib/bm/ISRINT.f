 
C*  8/30/91  K.Fujii       Calculate RSH distribution after ISR.
C*
      SUBROUTINE ISRINT(IPLOT,ROOTS)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8            ROOTS
C--
      PARAMETER       ( MXxDAT = 500 )
      COMMON   /RSHISR/ BMDATA(2,0:MXxDAT), BTE
      REAL*8            BMDATA            , BTE
C--
      REAL*8            RDATA(4,0:MXxDAT)
      DATA LOU /20/
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      xALF  = 1./137.
      xPI   = ACOS(-1.)
      AME   = 0.511E-3
      BTE   = (2*xALF/xPI)*(2*LOG(MAX(ROOTS,20*AME)/AME)-1)
C--
C  Loop over fraction.
C--
      NZ    = MXxDAT
      ZMN   = 1.D-6
      ZMX   = 1.D+0
      DZ    = (ZMX-ZMN)/NZ
C--
      PRB   = ( 1 + 3*BTE/4 )*ZMN
      DO 100 IZ = 0, NZ
         IF ( IZ.EQ.0 .OR. IZ.EQ.NZ ) THEN
            WGT = 1.D0/3
         ELSE IF ( MOD(IZ,2).EQ.0 ) THEN
            WGT = 2.D0/3
         ELSE
            WGT = 4.D0/3
         ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
         Z    = ZMN + IZ*DZ
         X    = Z**(1/BTE)
         RC   = BTE*(Z/X)*(1+3*BTE/4) - BTE*(1-X/2)
         DX   = DZ*(X/Z)/BTE
         XRS  = SQRT( 1 - X )
         ADD  = WGT*RC*DX
         PRB  = PRB + ADD
C--
         RDATA(1,IZ) = Z
         RDATA(2,IZ) = PRB
         RDATA(3,IZ) = XRS
         RDATA(4,IZ) = ADD/WGT/DX
100   CONTINUE
C--
C  Normalization.
C--
      DO 200 IZ = 0, NZ
         RDATA(2,IZ)  = RDATA(2,IZ)/PRB
         RDATA(4,IZ)  = RDATA(4,IZ)/PRB
         BMDATA(1,IZ) = RDATA(2,IZ)
         BMDATA(2,IZ) = RDATA(1,IZ)
200   CONTINUE
C--
C  Prepare topdraw file.
C--
      IF ( IPLOT.GT.0 ) THEN
         WRITE(LOU,*) 'NEW FRAME'
         WRITE(LOU,*) 'SET SCALE Y LOG'
         WRITE(LOU,*) 'SET LIMITS X 0 1'
         WRITE(LOU,*) 'SET ORDER X Y'
         DO 300 IZ = 0, NZ
            IF ( IPLOT.EQ.1 ) THEN
               WRITE(LOU,*) RDATA(3,IZ), RDATA(2,IZ)
            ELSE IF ( IPLOT.EQ.2 ) THEN
               WRITE(LOU,*) RDATA(3,IZ), RDATA(4,IZ)
            ENDIF
300      CONTINUE
         WRITE(LOU,*) 'JOIN'
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
