      SUBROUTINE MODGAM(PNAME,I,GAM)
 
      IMPLICIT    REAL*8 ( A-H, O-Z )
#include "../ss/include/SSPTAB.inc"
      CHARACTER*2 PNAME
      INTEGER  *4 I
      REAL     *8 GAM
C
C========< Entry Point >================================================
C
C--
C  Modify width.
C--
      IF ( PNAME.EQ.'SW' ) THEN
         IF ( I.GE.1 .AND. I.LE.2 ) THEN
            GMSW(I) = GAM
            PRINT *, ' MODGAM modified GMSW(',I,') to ', GAM 
         ELSE
							GO TO 9999
         ENDIF
      ELSE IF ( PNAME.EQ.'SZ' ) THEN
         IF ( I.GE.1 .AND. I.LE.4 ) THEN
            GMSZ(I) = GAM
            PRINT *, ' MODGAM modified GMSZ(',I,') to ', GAM 
         ELSE
							GO TO 9999
         ENDIF
      ELSE IF ( PNAME.EQ.'SF' ) THEN
         IF ( I.GE.1 .AND. I.LE.7 ) THEN
            GMSF(I) = GAM
            PRINT *, ' MODGAM modified GMSF(',I,') to ', GAM 
         ELSE
							GO TO 9999
         ENDIF
      ELSE IF ( PNAME.EQ.'SH' ) THEN
         IF ( I.GE.1 .AND. I.LE.4 ) THEN
            GMSH(I) = GAM
            PRINT *, ' MODGAM modified GMSH(',I,') to ', GAM 
         ELSE
							GO TO 9999
         ENDIF
      ELSE
							GO TO 9999
      ENDIF
C--
C  Normal end.
C--
      RETURN
C--
C  Error.
C--      
9999  PRINT *, ' >>>> MODGAM: ERROR '
      PRINT *, '    PNAME = ', PNAME, ' I = ', I, ' is invalid.'
      STOP
      END
