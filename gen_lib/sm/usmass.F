C***********************************************************************
C*
C*=====================================
C* Subroutine USMASS(MODE,PNAME,PMASS)
C*==============================-----==
C*
C* (Ppurpose)
C*    This subroutine allows users to put or get particle masses
C*    in or from /SMPTAB/ without directly referring to it.
C* (Inputs)
C*      MODE     : (R*4) : (1,2)=(get,put).
C*      PNAME    : (C*4) : particle name.
C* (Input/Output)
C*      PMASS    : (R*4) : particle mass.
C* (Relation)
C*    Calls no subroutines.
C* (Update Record)
C*   95/03/25  K.Fujii	Original version.
C*
C***********************************************************************
 
      SUBROUTINE USMASS(MODE,PNAME,PMASS)
 
      INTEGER  *4   MODE
      CHARACTER*4   PNAME
      REAL     *4   PMASS
C--
#include "../sm/include/SMPTAB.inc"
C
C========< Entry Point >================================================
C
C--
C  Decide the length of the name string.
C--
      L = INDEX(PNAME,' ') - 1 
C--
C  Branch on kinds of particles.
C--
      IF ( PNAME(:L).EQ.'E' ) THEN
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(1,2,1)
         ELSE
            AMF(1,2,1) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'MU' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(2,2,1)
         ELSE
            AMF(2,2,1) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'TAU' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(3,2,1)
         ELSE
            AMF(3,2,1) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'U' ) THEN
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(1,1,2)
         ELSE
            AMF(1,1,2) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'C' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(2,1,2)
         ELSE
            AMF(2,1,2) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'T' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(3,1,2)
         ELSE
            AMF(3,1,2) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'D' ) THEN
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(1,2,2)
         ELSE
            AMF(1,2,2) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'S' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(2,2,2)
         ELSE
            AMF(2,2,2) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'B' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMF(3,2,2)
         ELSE
            AMF(3,2,2) = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'W' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMW
         ELSE
            AMW = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'Z' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMZ
         ELSE
            AMZ = PMASS
         ENDIF
      ELSE IF ( PNAME(:L).EQ.'HSM' ) THEN 
         IF ( MODE.EQ.1 ) THEN
            PMASS = AMH
         ELSE
            AMH = PMASS
         ENDIF
      ELSE
         PRINT *, '>>> ERROR IN USMASS >>>'
         PRINT *, '  Particle ', PNAME(:L), ' is not supported.'
         PRINT *, '  Stop.'
         STOP
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
