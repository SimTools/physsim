CC**********************************************************************
C*
C*=========================================-----------------===
C* Subroutine UHQIJ3(AM,JPV,AMR,EPSL,QPR,Z,IPV,QIJ2,QIK2,WAT)
C*=========================================-----------------===
C*
C* (Purpose)
C*    Makes an appropriate choice of (ij) and (ik) for Dalitz's
C*    variables of 3-body phase space and return the generated
C*    values together with the Jacobian.
C* (Inputs)
C*       AM(i)          : (R*8) : mass of i
C*       JPV(1,r_1,r_2) : (I*4) : i    
C*          (2,r_1,r_2) : (I*4) : j   r_1 <--> (ij)
C*          (3,r_1,r_2) : (I*4) : k   r_2 <--> (ik)
C*       AMR(1,r)       : (R*8) : mass of resonance r
C*          (2,r)       : (R*8) : width of resonance r
C*          (3,r)       : (R*8) : m_(ij) + m_k  for r <--> (ij)
C*       EPSL           : (R*16) : width of gray area near poles
C*       QPR            : (R*16) : parent mass.
C*       Z(1-2)         : (R*16) : integration variables in (0,1)
C* (Output)
C*       IPV(1)   : (I*4) : i
C*          (2)   : (I*4) : j
C*          (3)   : (I*4) : k
C*       QIJ2     : (R*16) : (Q_ij)^2.
C*       QIK2     : (R*16) : (Q_ik)^2.
C*       WAT      : (R*16) : Jacobian weight.
C* (Relation)
C*    Invokes USORTR in FORTLIB and UHQLIM and UHQ2BW in ttz_lib.
C* (Update Record)
C*   95/04/07  K.Fujii	Original version using BW parametrization.
C*
CC**********************************************************************

      SUBROUTINE UHQIJ3(AM,JPV,AMR,EPSL,QPR,Z,IPV,QIJ2,QIK2,WAT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      INTEGER*4  JPV(3,3,3), IPV(3)
      REAL*8  AM(*), AMR(3,3)
      REAL*8  EPSL, QPR, Z(2), QIJ2, QIK2, WAT
C--
C  Local variables.
C--
      INTEGER*4  ISHF(3)
C--
      DATA NCALL /  0 /
C
C========< Entry Point >================================================
C
C--
C  Reset event weight.
C--
      WAT = 1
C--
C  Sort possible resonances.
C     AMR(1,r') = AM_ij
C        (2,r') = GM_ij
C        (3,r') = AM_ij + AM_k
C  where r = 1, 2, 3 correspond to 3 combinations of
C  daughter particles (i,j,k). Notice that i, j, and k point
C  to AM and PV and not necessarily congiguous.
C     ISHF(r) = r' for the r-th resonance from the lightest.
C--
      CALL USORTR(3,3,3,AMR(1,1),ISHF(1))
C--
C  Check if there are more than 2 poles.
C--
      IF ( AMR(3,ISHF(3)).GT.QPR+EPSL ) THEN
C--
C  -----------------
C  Less than 3 poles 
C  -----------------
C  Choose the pairs corresponding to the lightest two as 
C  independent variables.
C     JPV(1,r'1,r'2) = i 
C        (2,r'1,r'2) = j
C        (3,r'1,r'2) = k
C  where r'1 is a (ij) resonance and r'2 is a (ik) resonance. 
C--
         IPV(1) = JPV(1,ISHF(1),ISHF(2))   
         IPV(2) = JPV(2,ISHF(1),ISHF(2))
         IPV(3) = JPV(3,ISHF(1),ISHF(2))
C--
         AM_A = AM(IPV(1))
         AM_B = AM(IPV(2))
         AM_C = AM(IPV(3))
C-- ab.
         QIJMN = AM_A + AM_B
         QIJMX = QPR - AM_C
         CALL UHQ2BW(DBLE(AMR(1,ISHF(1))),DBLE(AMR(2,ISHF(1))),
     .                                QIJMN,QIJMX,Z(1),QIJ2,WT)
         QIJ   = SQRT(QIJ2)
         WAT   = WAT*WT
C-- ac.
         CALL UHQLIM(QPR,AM_A,AM_B,AM_C,QIJ,QIKMN,QIKMX)
         CALL UHQ2BW(DBLE(AMR(1,ISHF(2))),DBLE(AMR(2,ISHF(2))),
     .                                QIKMN,QIKMX,Z(2),QIK2,WT)
         WAT   = WAT*WT
C--
      ELSE
C--
C  -------
C  3 poles 
C  -------
C  Always choose the pair giving the lightest (bc) as an in dependent
C  variables, while the other should be selected from (ab) or (ac), 
C  depending on the region in the phase space.
C--
         XR  = AMR(1,ISHF(2)) - EPSL
         YR  = AMR(1,ISHF(3)) - EPSL
C--
         IA  = JPV(1,ISHF(2),ISHF(3))
         IB  = JPV(2,ISHF(2),ISHF(3))
         IC  = JPV(3,ISHF(2),ISHF(3))
C--
         AM_A = AM(IA)
         AM_B = AM(IB)
         AM_C = AM(IC)
C--
C  Do not allow the double resonance of the heaviest two.
C  This is a restriction to the parameter space for simplicity of
C  the program (and is due to my laziness).
C--         
         CALL UHQLIM(QPR,AM_A,AM_B,AM_C,XR,YRMN,YRMX)
         CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YR,XRMN,XRMX)
C--         
         IF ( YR.LE.YRMX .OR. XR.LE.XRMX ) THEN
            PRINT *, '>>>> GETXPH: CAUTION '
            PRINT *, '  Kinematics inappropriate for this set of ',
     .                 'parameters.'
            PRINT *, '     QPR = ', QPR
            PRINT *, '     AM  = ', AM_A, AM_B, AM_C
            PRINT *, '     AMR = ', AMR(1,ISHF(1)),
     .                              AMR(1,ISHF(2)),
     .                              AMR(1,ISHF(3))
            PRINT *, '  Ignore this point.'
            PRINT *, '  When you got this message too often, ',
     .                 'the result might be wrong.'
            WAT = 0
            RETURN
         ENDIF
C--
C  The point accepted.
C--
         XMN   = AM_A + AM_B
         XMX   = QPR - AM_C         
         XCT   = (XR+XRMX)/2
C--         
C  Branch on kinematical regions.
C--
         IF ( XMN+Z(1)*(XMX-XMN).LT.XCT ) THEN
C--
C  Region I: choose (ac) and (bc) as independent.
C--
            Z1  = Z(1)*(XMX-XMN)/(XCT-XMN)
            WAT = WAT*(XMX-XMN)/(XCT-XMN)
C--
            IPV(1) = IC  
            IPV(2) = IA
            IPV(3) = IB
C--            
            YMN = AM_A + AM_C
            YMX = QPR - AM_B
            CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YMN,XCMN,XCMX)
            XC  = XCMN
            YC  = YMN
C-- ac.
            IF ( XC.LE.XCT ) THEN
               QIJMN = YMN
            ELSE
               CALL UHQLIM(QPR,AM_A,AM_B,AM_C,XCT,YCTMN,YCTMX)
               QIJMN = YCTMN
            ENDIF
            QIJMX = YMX
            CALL UHQ2BW(DBLE(AMR(1,ISHF(3))),DBLE(AMR(2,ISHF(3))),
     .                                     QIJMN,QIJMX,Z1,QIJ2,WT)
            QIJ   = SQRT(QIJ2)
            WAT   = WAT*WT
C-- bc.
            CALL UHQLIM(QPR,AM_C,AM_A,AM_B,QIJ,QIKMN,QIKMX)
            IF ( QIJ.GT.YCTMN .AND. QIJ.LT.YCTMX ) THEN
               QIKMN = AM_A**2 + AM_B**2 + AM_C**2
     .                + (QPR-XCT)*(QPR+XCT) - QIJ**2
               QIKMN = SQRT(MAX(QIKMN,0.D0))                
            ENDIF
            CALL UHQ2BW(DBLE(AMR(1,ISHF(1))),DBLE(AMR(2,ISHF(1))),
     .                                   QIKMN,QIKMX,Z(2),QIK2,WT)
            WAT   = WAT*WT
C--
         ELSE
C--
C  Region II: choose (ab) and (bc) as independent.
C--
            Z1  = ((XMX-XMN)*Z(1)-(XCT-XMN))/(XMX-XCT)
            WAT = WAT*(XMX-XMN)/(XMX-XCT)
C--
            IPV(1) = IB
            IPV(2) = IA
            IPV(3) = IC
C-- ab.
            QIJMN  = XCT
            QIJMX  = XMX
            CALL UHQ2BW(DBLE(AMR(1,ISHF(2))),DBLE(AMR(2,ISHF(2))),
     .                                     QIJMN,QIJMX,Z1,QIJ2,WT)
            QIJ   = SQRT(QIJ2)
            WAT   = WAT*WT
C-- bc.
            CALL UHQLIM(QPR,AM_B,AM_A,AM_C,QIJ,QIKMN,QIKMX)
            CALL UHQ2BW(DBLE(AMR(1,ISHF(1))),DBLE(AMR(2,ISHF(1))),
     .                                   QIKMN,QIKMX,Z(2),QIK2,WT)
            WAT   = WAT*WT
         ENDIF
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
