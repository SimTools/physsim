CC********************************************************************CC
C*                                                                    *C
C*==========================================================          *C
C* Subroutine UCONTR(FNCT,NCT,CONTUR,XYLIM,XY0,DXY0,LOUNIT)           *C
C*==========================================================          *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*     Produce TOPDRAW data for a contour plot. FNCT is a             *C
C*     external function of R*4 variables (X,Y).                      *C
C* (Inputs)                                                           *C
C*     FNCT       : (R*4) ; a R*4 function of R*4 (X,Y).              *C
C*     NCT        : (I*4) ; # contours to draw.                       *C
C*     CONTUR(I)  : (R*4) ; FNCT value corresponding to I-th contour. *C
C*     XYLIM(1,1) : (R*4) ; X-min --+                                 *C
C*          (2,1) : (R*4) ; X-max   +-- of the window.                *C
C*          (1,2) : (R*4) ; Y-min   !                                 *C
C*          (2,2) : (R*4) ; Y-max --+                                 *C
C*     XY0(1)     : (R*4) ; X --+--- of a point in the window.        *C
C*        (2)     : (R*4) ; Y --+                                     *C
C*     DXY0(1)    : (R*4) ; X-error.                                  *C
C*         (2)    : (R*4) ; Y-error.                                  *C
C*     LOUNIT     : (I*4) ; logical unit number for TDR data.         *C
C* (Output)                                                           *C
C*    TOPDRAW data is output to LOUNIT.                               *C
C* (Relation)                                                         *C
C*    Calls UGRADX.                                                   *C
C* (Update Record)                                                    *C
C*    02/25/90  K.Fujii            Original version.                  *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UCONTR(FNCT,NCT,CONTUR,XYLIM,XY0,DXY0,LOUNIT)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      INTEGER*4  NCT
      REAL*8     CONTUR(NCT), XYLIM(2,2), XY0(2), DXY0(2)
      REAL*8     SKL(2), XY(2), DXY(2), DF(2)
C
C     EPS for controlling maximum scaled step size.
C     EPSL for convergence check.
C     DFIMN is the minimum d(phi)/d(step).
C     DFIMX is the maximum d(phi)/d(step).
C     STEPMN is the minimum scaled step length.
C
      DATA       EPS    / 5.E-1 /   EPSL  / 5.E-5 /
      DATA       DFIMN  / 1.E-1 /   DFIMX / 5.E-1 /
     .           STEPMN / 2.E-2 /
      DATA       NITCUT / 10 /
      DATA       MXxPNT / 2000 /
      DATA       NCALL  / 0    /
      EXTERNAL   FNCT
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         xPI   = ACOS(-1.E0)
         x2PI  = 2*xPI
      ENDIF
C--
C  Define window of the contour plot.
C--
      XMN    = XYLIM(1,1)
      XMX    = XYLIM(2,1)
      YMN    = XYLIM(1,2)
      YMX    = XYLIM(2,2)
      XSKL   = ABS(XYLIM(2,1)-XYLIM(1,1))
      YSKL   = ABS(XYLIM(2,2)-XYLIM(1,2))
      SKL(1) = XSKL
      SKL(2) = YSKL
C>>> 91/06/01
      STEPMX = SQRT((DXY0(1)/XSKL)**2 + (DXY0(2)/YSKL)**2)
C>>>
C--
C  Start find contours.
C--
      WRITE(LOUNIT,*) 'SET LIMIT X ', XMN, XMX
      WRITE(LOUNIT,*) 'SET LIMIT Y ', YMN, YMX
C--
      DO 5000 ICT = 1, NCT
         NTRTOT = 0
C--
         WRITE(LOUNIT,*) '( ICT =', ICT
         WRITE(LOUNIT,*) 'SET ORDER X Y'
         PRINT *, ' Start drawing ', ICT, '-th contour.'
C--
C  Decide starting point.
C--
         FCT    = CONTUR(ICT)
         XY(1)  = XY0(1)
         XY(2)  = XY0(2)
         DXY(1) = DXY0(1)
         DXY(2) = DXY0(2)
         F0     = FNCT(XY(1),XY(2))
C>>> 91/05/30
C        FTEST  = EPSL*ABS(FCT-F0)
         FTEST  = MAX(EPSL*ABS(FCT-F0),EPSL)
C>>>
         NIT0   = 0
1        CALL UGRADX(FNCT,FCT,FTEST,SKL,XY,DXY,DF,NTRY)
         NTRTOT = NTRTOT + IABS(NTRY)
C--
         IF ( NTRY.LT.0 .OR. XY(1).LT.XMN .OR. XY(1).GT.XMX .OR.
     .        XY(2).LT.YMN .OR. XY(2).GT.YMX ) THEN
            NIT0   = NIT0 + 1
            XY(1)  = XY0(1) + DXY(1)
            XY(2)  = XY0(2) + DXY(2)
            DXY(1) = DXY(1)/2
            DXY(2) = DXY(2)/2
            IF ( NIT0.LT.10 )                    GO TO 1
            PRINT *, ' Failed to get starting point.'
            PRINT *, ' STOP in CONT.'
C           STOP
                                                 GO TO 5000
         ENDIF
         PRINT *, ' NTRY = ', NTRY, ' to get starting point.'
C--
C  Starting point decided. Now decide step size.
C--
         X0     = XY(1)
         Y0     = XY(2)
         DFX    = DF(1)*XSKL
         DFY    = DF(2)*YSKL
         FI0    = MOD(ATAN2(DFY,DFX)+4*x2PI,x2PI)
         STEP0  = EPS*SQRT(((X0-XY0(1))/XSKL)**2
     .                   + ((Y0-XY0(2))/YSKL)**2 )
C>>> 91/05/30
         DFI    = 0
         STEP0  = MAX(STEP0,STEPMN)
C>>>
C--
C  Now trace this contour.
C--
         ITRY  = 0
2        ITRY  = ITRY + 1
C--
C  Recall starting point.
C--
         XY(1) = X0
         XY(2) = Y0
         FI    = FI0
C--
         WRITE(LOUNIT,*) XY(1), XY(2)
C--
C  Looping starts here.
C--
         STEP   = STEP0
         NPNT   = 0
         DELFI  = 0
         NIT1   = 0
         NIT2   = 0
C--
3        NPNT   = NPNT + 1
C--
         XS     = XY(1)
         YS     = XY(2)
         DFXS   = DFX
         DFYS   = DFY
         FIS    = FI
C--
         ADF    = SQRT(DFX**2+DFY**2)
C>>> 91/06/01
C        IF ( ITRY.EQ.1 ) THEN
C          DELX   =-XSKL*STEP*DFY/ADF
C          DELY   = YSKL*STEP*DFX/ADF
C        ELSE
C          DELX   = XSKL*STEP*DFY/ADF
C          DELY   =-YSKL*STEP*DFX/ADF
C        ENDIF
C
         IF ( ITRY.EQ.1 ) THEN
           DELX   =-STEP*DFY/ADF
           DELY   = STEP*DFX/ADF
         ELSE
           DELX   = STEP*DFY/ADF
           DELY   =-STEP*DFX/ADF
         ENDIF
         DELX = DELX - DFIDS*STEP/2*DELY
         DELY = DELY + DFIDS*STEP/2*DELX
         DELX = XSKL*DELX
         DELY = YSKL*DELY
C>>>
         XY(1)  = XY(1) + DELX
         XY(2)  = XY(2) + DELY
C>>> 91/05/30
C        PRINT '(F12.5,4F12.3)', STEP, XY(1), XY(2), DELX, DELY
C>>>
         DXY(1) = DELX
         DXY(2) = DELY
         CALL UGRADX(FNCT,FCT,FTEST,SKL,XY,DXY,DF,NTRY)
         NTRTOT = NTRTOT + IABS(NTRY)
         IF ( IABS(NTRY).GT.10 ) THEN
CX          PRINT *, ' NPNT = ', NPNT, ' NTRY = ', NTRY
         ENDIF
C--
C  Check if this step is reasonalble.
C--
         DFX    = XSKL*DF(1)
         DFY    = YSKL*DF(2)
C>>> 91/05/30
CX       IF ( DFX.LE.1.D-50 .AND. DFY.LE.1.D-50 ) THEN
CX          DFX = DFXS
CX          DFY = DFYS
CX       ENDIF
         ADF    = SQRT(DFX**2+DFY**2)
         IF ( ADF.GT.1.D-50 ) THEN
            FI     = MOD(ATAN2(DFY,DFX)+4*x2PI,x2PI)
         ENDIF
C>>>
         DFI    = MOD(FI-FIS+4*x2PI,x2PI)
         IF ( DFI.GT.xPI ) DFI = DFI - x2PI
C>>> 91/05/30
C        ADFI   = ABS(DFI)
         ADFI   = MAX(ABS(DFI),1.D-50)
         DFIDS  = DFI/STEP
C>>>
C        PRINT *, ' ADFI NTRY NIT1 NIT2 = ', ADFI, NTRY, NIT1, NIT2
C>>>
         IF ( ADFI.GT.DFIMX .OR. NTRY.LT.0 ) THEN
            NIT1  = NIT1 + 1
            IF ( NIT1.GT.NITCUT ) THEN
               IF ( NTRY.LT.0 ) THEN
                  PRINT *, ' BLACK HOLE at NPNT = ', NPNT, ' in CONT.'
                  PRINT *, '    NIT1, NIT2 = ', NIT1, NIT2
                  PRINT *, ' Give up continuing and STOP.'
C>>> 91/05/30
C                 STOP
                  WRITE(LOUNIT,*) 'PLOT'
                  WRITE(LOUNIT,*) 'JOIN'
                  IF ( ITRY.LT.2 )               GO TO 2
                  RETURN
C>>>
               ELSE
                                                 GO TO 4
               ENDIF
            ENDIF
            XY(1) = XS
            XY(2) = YS
            DFX   = DFXS
            DFY   = DFYS
            FI    = FIS
            NPNT  = NPNT - 1
            STEP  = MAX(STEP*MIN(0.3,DFIMX/ADFI),STEPMN)
C>>> 91/06/01
C           PRINT *, ' Reduce STEP = ', STEP
C>>>
                                                 GO TO 3
         ELSE IF ( ADFI.LT.DFIMN .AND. NIT2.LT.NITCUT .AND.
     .             ABS(STEP-STEP0).GT.1.D-20 ) THEN
            NIT2  = NIT2 + 1
            XY(1) = XS
            XY(2) = YS
            DFX   = DFXS
            DFY   = DFYS
            FI    = FIS
            NPNT  = NPNT - 1
C>>> 91/06/01
C           STEP  = MIN(STEP*MAX(3.,DFIMN/ADFI),STEP0)
            STEP  = MIN(STEP*MAX(3.,DFIMN/ADFI),STEPMX)
C           PRINT *, ' Increase STEP = ', STEP
C>>>
                                                 GO TO 3
         ENDIF
C--
C  Point accepted.
C--
4        NIT1 = 0
         NIT2 = 0
C>>> 91/05/30
C        PRINT '(F12.5,4F12.3)', STEP, XY(1), XY(2), XS, YS
C        WRITE(LOUNIT,*) XY(1), XY(2)
         EXMX = ABS(XYLIM(2,1)-XYLIM(1,1))/10
         EYMX = ABS(XYLIM(2,2)-XYLIM(1,2))/10
         IF ( XY(1).GE.XMN-EXMX .AND. XY(1).LE.XMX+EXMX .AND.
     .        XY(2).GE.YMN-EYMX .AND. XY(2).LE.YMX+EYMX .AND.
     .        NTRY.GT.0 ) THEN
            WRITE(LOUNIT,*) XY(1), XY(2)
         ENDIF
C>>>
C--
C  Boundary check.
C--
         IF ( XY(1).LT.XMN .OR. XY(1).GT.XMX .OR.
     .        XY(2).LT.YMN .OR. XY(2).GT.YMX ) THEN
            IF ( ITRY.GE.2 )                     GO TO 40
            WRITE(LOUNIT,*) 'PLOT'
            WRITE(LOUNIT,*) 'JOIN'
                                                 GO TO 2
         ENDIF
C--
C  Check if this loop is completed.
C--
         DELFI  = DELFI + DFI
         IF ( ABS(DELFI).GT.x2PI + 1.E-4 )       GO TO 40
C        DIST   = SQRT(((XY(1)-X0)/XSKL)**2+((XY(2)-Y0)/YSKL)**2)
C        IF ( ABS(DELFI).GT.x2PI + 1.E-4 .AND.
C    .        DIST.LT.STEP0 )                    GO TO 40
C--
C  If not, try next point.
C--
         IF ( NPNT.LT.MXxPNT )                   GO TO 3
         PRINT *, ' Step size too small in CONT.'
         PRINT *, ' NPNT = ', NPNT, ' > MXxPNT = ', MXxPNT
C--
C  Contour completed.
C--
40       WRITE(LOUNIT,*) 'PLOT'
         WRITE(LOUNIT,*) 'JOIN'
         PRINT *, ICT, '-th contour completed. NTRTOT = ', NTRTOT
5000  CONTINUE
C--
C  That's it.
C--
      RETURN
      END
