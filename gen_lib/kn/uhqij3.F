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
C*       AM(i)          : (R*4) : mass of i
C*       JPV(1,r_1,r_2) : (I*4) : i    
C*          (2,r_1,r_2) : (I*4) : j   r_1 <--> (ij)
C*          (3,r_1,r_2) : (I*4) : k   r_2 <--> (ik)
C*       AMR(1,r)       : (R*4) : mass of resonance r
C*          (2,r)       : (R*4) : width of resonance r
C*          (3,r)       : (R*4) : m_(ij) + m_k  for r <--> (ij)
C*       EPSL           : (R*8) : width of gray area near poles
C*       QPR            : (R*8) : parent mass.
C*       Z(1-2)         : (R*8) : integration variables in (0,1)
C* (Output)
C*       IPV(1)   : (I*4) : i
C*          (2)   : (I*4) : j
C*          (3)   : (I*4) : k
C*       QIJ2     : (R*8) : (Q_ij)^2.
C*       QIK2     : (R*8) : (Q_ik)^2.
C*       WAT      : (R*8) : Jacobian weight.
C* (Relation)
C*    Invokes USORTR in FORTLIB and UHQLIM and UHQ2BW in ttz_lib.
C* (Update Record)
C*   95/04/07  K.Fujii	Original version using BW parametrization.
C*   03/09/16  K.Fujii  Modified to handle general 3-pole cases.
C*
CC**********************************************************************

      SUBROUTINE UHQIJ3(AM,JPV,AMR,EPSL,QPR,Z,IPV,QIJ2,QIK2,WAT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
C--
C  Dummy arguments.
C--
      INTEGER*4  JPV(3,3,3), IPV(3)
      REAL   *4  AM(*), AMR(3,3)
      REAL   *8  EPSL, QPR, Z(2), QIJ2, QIK2, WAT
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
      IRG = 0
      IF ( AMR(3,ISHF(3)).GT.QPR+MIN(EPSL,50*AMR(2,ISHF(3))) ) THEN
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
C  Always choose the pair giving the lightest (bc) as an independent
C  variables, while the other should be selected from (ab) or (ac), 
C  depending on the region in the phase space.
C--
         XR  = AMR(1,ISHF(2))
         YR  = AMR(1,ISHF(3)) - MIN(EPSL,10*AMR(2,ISHF(3)))
C--
         IA  = JPV(1,ISHF(2),ISHF(3))
         IB  = JPV(2,ISHF(2),ISHF(3))
         IC  = JPV(3,ISHF(2),ISHF(3))
C--
         AM_A = AM(IA)
         AM_B = AM(IB)
         AM_C = AM(IC)
C--
         XMN = AM_A + AM_B
         XMX = QPR  - AM_C         
C--
C  Check if the heaviest two may simultaneously be resonant.
C--         
         CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YR,XRMN,XRMX)
C--        
         IF ( XR.GE.XRMN .AND. XR.LE.XRMX ) THEN
            IRG = IRG + 10
C--
C  The heaviest two may simultaneously be on-shell.
C--
            QIJR = AMR(1,ISHF(1))
            XP = AM_A**2 + AM_B**2 + AM_C**2
     .           + (QPR-YR)*(QPR+YR) - QIJR**2
            XP = SQRT(MAX(XP,0.D0))
            YP = YR
C--
            XQ = XR
            YQ = AM_A**2 + AM_B**2 + AM_C**2
     .           + (QPR-XR)*(QPR+XR) - QIJR**2
            YQ = SQRT(MAX(YQ,0.D0))
C--
            XCT = (XR+XP)/2
            IF ( XP.LT.XR ) THEN
               ICASE = +1
            ELSE
               ICASE = -1
            ENDIF
            YCT = (YR + (YR+YQ)/2)/2
         ELSE
C--
C  The heaviest two cannot be simultaneously on-shell.
C--
            IF ( XR.GT.XRMX ) THEN
               XCT   = (XR+XRMX)/2
               ICASE = +1
            ELSE
               XCT   = (XR+XRMN)/2
               ICASE = -1
            ENDIF
            YCT = -9999
         ENDIF
C--         
C  Branch on kinematical regions.
C--
         IF ( (XMN+Z(1)*(XMX-XMN))*ICASE.LT.XCT*ICASE ) THEN
C--
C  Region I: choose (ac) and (bc) as independent.
C--
            IRG = IRG + 1
C--
            IF ( ICASE.GT.0 ) THEN
               Z1  = Z(1)*(XMX-XMN)/(XCT-XMN)
               WAT = WAT*(XMX-XMN)/(XCT-XMN)
            ELSE
               Z1  = ((XMX-XMN)*Z(1)-(XCT-XMN))/(XMX-XCT)
               WAT = WAT*(XMX-XMN)/(XMX-XCT)
            ENDIF
C--
            YMN = AM_A + AM_C
            YMX = QPR  - AM_B
            IPV(1) = IC  
            IPV(2) = IA
            IPV(3) = IB
C--            
            CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YMN,XCMN,XCMX)
            XC  = XCMN
            YC  = YMN
            CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YMX,XDMN,XDMX)
            XD  = XDMN
            YD  = YMX
            CALL UHQLIM(QPR,AM_A,AM_B,AM_C,XCT,YCTMN,YCTMX)
C-- ac.
            IF ( XC*ICASE.LE.XCT*ICASE ) THEN
               QIJMN = YMN
            ELSE
               QIJMN = YCTMN
            ENDIF
C--
            IF ( XD*ICASE.LE.XCT*ICASE ) THEN
               QIJMX = YMX
            ELSE
               QIJMX = YCTMX
            ENDIF
C--
            CALL UHQ2BW(DBLE(AMR(1,ISHF(3))),DBLE(AMR(2,ISHF(3))),
     .                                     QIJMN,QIJMX,Z1,QIJ2,WT)
            QIJ   = SQRT(QIJ2)
            WAT   = WAT*WT
C-- bc.
            CALL UHQLIM(QPR,AM_C,AM_A,AM_B,QIJ,QIKMN,QIKMX)
C--
            IF ( QIJ.GT.YCTMN .AND. QIJ.LT.YCTMX ) THEN
               QIKCT = AM_A**2 + AM_B**2 + AM_C**2
     .                + (QPR-XCT)*(QPR+XCT) - QIJ**2
               QIKCT = SQRT(MAX(QIKCT,0.D0))                
               IF ( ICASE.LT.0 ) THEN
                  QIKMX = QIKCT
               ELSE
                  QIKMN = QIKCT
               ENDIF 
            ENDIF
C--
            CALL UHQ2BW(DBLE(AMR(1,ISHF(1))),DBLE(AMR(2,ISHF(1))),
     .                                   QIKMN,QIKMX,Z(2),QIK2,WT)
            WAT   = WAT*WT
C--
         ELSE
C--
C  Region II or III.
C--
            IF ( ICASE.GT.0 ) THEN
               Z1  = ((XMX-XMN)*Z(1)-(XCT-XMN))/(XMX-XCT)
               WAT = WAT*(XMX-XMN)/(XMX-XCT)
            ELSE
               Z1  = Z(1)*(XMX-XMN)/(XCT-XMN)
               WAT = WAT*(XMX-XMN)/(XCT-XMN)
            ENDIF
C--
            IF ( YCT.LT.0.D0 ) THEN     
C--
C  Region II (No Region III): choose (ab) and (bc) as independent.
C--
               IRG = IRG + 2          
C--
               IPV(1) = IB
               IPV(2) = IA
               IPV(3) = IC
C-- ab.
               IF ( ICASE.GT.0 ) THEN
                  QIJMN  = XCT
                  QIJMX  = XMX
               ELSE
                  QIJMN  = XMN
                  QIJMX  = XCT
               ENDIF
               CALL UHQ2BW(DBLE(AMR(1,ISHF(2))),DBLE(AMR(2,ISHF(2))),
     .                                        QIJMN,QIJMX,Z1,QIJ2,WT)
               QIJ   = SQRT(QIJ2)
               WAT   = WAT*WT
C-- bc.
               CALL UHQLIM(QPR,AM_B,AM_A,AM_C,QIJ,QIKMN,QIKMX)
               CALL UHQ2BW(DBLE(AMR(1,ISHF(1))),DBLE(AMR(2,ISHF(1))),
     .                                   QIKMN,QIKMX,Z(2),QIK2,WT)
               WAT   = WAT*WT
            ELSE
C--
C  Region II or III.
C--
               YMN = AM_A + AM_C
               YMX = QPR  - AM_B
C--            
               CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YMN,XCMN,XCMX)
               XC  = XCMN
               YC  = YMN
               CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YMX,XDMN,XDMX)
               XD  = XDMN
               YD  = YMX
               CALL UHQLIM(QPR,AM_A,AM_B,AM_C,XMN,YEMN,YEMX)
               XE  = XMN
               YE  = YEMN
               CALL UHQLIM(QPR,AM_A,AM_B,AM_C,XMX,YFMN,YFMX)
               XF  = XMX
               YF  = YFMX
C--
               CALL UHQLIM(QPR,AM_A,AM_C,AM_B,YCT,XCTMN,XCTMX)
               CALL UHQLIM(QPR,AM_A,AM_B,AM_C,XCT,YCTMN,YCTMX)
C--
               IF ( XC*ICASE.GT.XCT*ICASE ) THEN
                  QIJMN = YMN
               ELSE
                  QIJMN = YCTMN
               ENDIF
C--
               IF ( XD*ICASE.GT.XCT*ICASE ) THEN
                  QIJMX = YMX
               ELSE
                  QIJMX = YCTMX
               ENDIF
C--
               IF ( (QIJMN+Z1*(QIJMX-QIJMN))*ICASE.LT.YCT*ICASE ) THEN
C--
C  Region II: choose (ab) and (bc) as independent.
C--
                  IRG = IRG + 2
C--          
                  IF ( ICASE.GT.0 ) THEN
                     Z1  = Z1*(QIJMX-QIJMN)/(YCT-QIJMN)
                     WAT = WAT*(QIJMX-QIJMN)/(YCT-QIJMN)
                  ELSE
                     Z1  = ((QIJMX-QIJMN)*Z1-(YCT-QIJMN))/(QIJMX-YCT)
                     WAT = WAT*(QIJMX-QIJMN)/(QIJMX-YCT)
                  ENDIF
C--
                  IPV(1) = IB
                  IPV(2) = IA
                  IPV(3) = IC
C-- ab.
                  IF ( ICASE.GT.0 ) THEN
                     QIJMN  = XCT
                     IF ( YF.LT.YCT ) THEN
                        QIJMX = XMX
                     ELSE
                        QIJMX = XCTMX
                     ENDIF
                  ELSE
                     IF ( YE.GT.YCT ) THEN
                        QIJMN  = XMN
                     ELSE
                        QIJMN  = XCTMN
                     ENDIF
                     QIJMX  = XCT
                  ENDIF
C--
                  CALL UHQ2BW(DBLE(AMR(1,ISHF(2))),DBLE(AMR(2,ISHF(2))),
     .                                           QIJMN,QIJMX,Z1,QIJ2,WT)
                  QIJ   = SQRT(QIJ2)
                  WAT   = WAT*WT
C-- bc.
                  CALL UHQLIM(QPR,AM_B,AM_A,AM_C,QIJ,QIKMN,QIKMX)

                  QIKCT = AM_A**2 + AM_B**2 + AM_C**2
     .                   + (QPR-YCT)*(QPR+YCT) - QIJ**2
                  QIKCT = SQRT(MAX(QIKCT,0.D0))
C--
                  IF ( ICASE.GT.0.D0 ) THEN
                     QIKMN = QIKCT
                  ELSE
                     QIKMX = QIKCT
                  ENDIF
C--
                  CALL UHQ2BW(DBLE(AMR(1,ISHF(1))),DBLE(AMR(2,ISHF(1))),
     .                                      QIKMN,QIKMX,Z(2),QIK2,WT)
                  WAT   = WAT*WT
               ELSE
C--
C  Region III: choose (ab) and (ac) as independent.
C--
                  IRG = IRG + 3         
C--
                  IF ( ICASE.GT.0 ) THEN
                     Z1  = ((QIJMX-QIJMN)*Z1-(YCT-QIJMN))/(QIJMX-YCT)
                     WAT = WAT*(QIJMX-QIJMN)/(QIJMX-YCT)
                  ELSE
                     Z1  = Z1*(QIJMX-QIJMN)/(YCT-QIJMN)
                     WAT = WAT*(QIJMX-QIJMN)/(YCT-QIJMN)
                  ENDIF
C--
                  IPV(1) = IA
                  IPV(2) = IB
                  IPV(3) = IC
C-- ab.
                  IF ( ICASE.GT.0 ) THEN
                     QIJMN  = XCT
                     IF ( YF.GE.YCT ) THEN
                        QIJMX = XMX
                     ELSE
                        QIJMX = XCTMX
                     ENDIF
                  ELSE
                     IF ( YE.LE.YCT ) THEN
                        QIJMN  = XMN
                     ELSE
                        QIJMN  = XCTMN
                     ENDIF
                     QIJMX  = XCT
                  ENDIF
C--
                  CALL UHQ2BW(DBLE(AMR(1,ISHF(2))),DBLE(AMR(2,ISHF(2))),
     .                                           QIJMN,QIJMX,Z1,QIJ2,WT)
                  QIJ   = SQRT(QIJ2)
                  WAT   = WAT*WT
C-- ac.
                  CALL UHQLIM(QPR,AM_A,AM_B,AM_C,QIJ,QIKMN,QIKMX)
C--
                  IF ( ICASE.GT.0.D0 ) THEN
                     QIKMN = YCT
                  ELSE
                     QIKMX = YCT
                  ENDIF
C--
                  CALL UHQ2BW(DBLE(AMR(1,ISHF(3))),DBLE(AMR(2,ISHF(3))),
     .                                         QIKMN,QIKMX,Z(2),QIK2,WT)
                  WAT   = WAT*WT
C--
               ENDIF
            ENDIF
         ENDIF
      ENDIF
#ifdef __DEBUG__
      IF ( IRG.GT.10 ) THEN
            PRINT *, '>>>> GETXPH:  '
            PRINT *, '     XR   = ', XR
            PRINT *, '     XRMN = ', XRMN
            PRINT *, '     XRMX = ', XRMX
            PRINT *, '     QPR = ', QPR
            PRINT *, '     AM  = ', AM_A, AM_B, AM_C
            PRINT *, '     AMR = ', AMR(1,ISHF(1)),
     .                              AMR(1,ISHF(2)),
     .                              AMR(1,ISHF(3))
            PRINT *, '     IRG = ', IRG
            PRINT *, '     IPV = ', IPV
            PRINT *, '     QIJ = ', SQRT(QIJ2)
            PRINT *, '     QIK = ', SQRT(QIK2)
      ENDIF
#endif
C--
C  That's it.
C--
      RETURN
      END
