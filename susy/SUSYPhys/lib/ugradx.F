 
CC********************************************************************CC
C*                                                                    *C
C*=======================================================             *C
C* Subroutine UGRADX(FNCT,FCT,FTEST,SKL,XY,DXY,DF,NTRY)               *C
C*=======================================================             *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    UGRADX looks for (X,Y) at which ABS(FNCT(X,Y)-FCT)<FTEST,       *C
C*    moving along the gradient vector from (XY(1),XY(2)).            *C
C*    (DXY(1),DXY(2)) gives initial step size.                        *C
C*    When such (X,Y) is found, the coordinates are saved in          *C
C*    XY, the last step in (DXY), and the gradient vector at          *C
C*    (X,Y) in DF.                                                    *C
C* (Inputs)                                                           *C
C*     FNCT       : (R*4) ; a R*4 function of R*4 (X,Y).              *C
C*     FCT        : (I*4) ; function value on the contour.            *C
C*     FTEST      : (R*4) ; TOLERANCE ON THE FUNCTION VALUE.          *C
C*     SKL(1)     : (R*4) ; X-scale.                                  *C
C*        (2)     : (R*4) ; Y-scale.                                  *C
C*     XY (1)     : (R*4) ; X of the current position.                *C
C*        (2)     : (R*4) ; Y of the current position.                *C
C*     DXY(1)     : (R*4) ; X step.                                   *C
C*        (2)     : (R*4) ; Y step.                                   *C
C* (Outputs)                                                          *C
C*     XY (1)     : (R*4) ; X of the next position.                   *C
C*        (2)     : (R*4) ; Y of the next position.                   *C
C*     DXY(1)     : (R*4) ; accepted X step.                          *C
C*        (2)     : (R*4) ; accepted Y step.                          *C
C*     DF (1)     : (R*4) ; @F/@X at the next position.               *C
C*        (2)     : (R*4) ; @F/@Y at the next position.               *C
C* (Relation)                                                         *C
C*    Calls no subroutines.                                           *C
C* (Update Record)                                                    *C
C*    02/25/90   K.Fujii                 Original version.            *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UGRADX(FNCT,FCT,FTEST,SKL,XY,DXY,DF,NTRY)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      INTEGER*4  NTRY
      REAL*8     FCT, FTEST, SKL(2), XY(2), DXY(2), DF(2)
C
C     EPS for derivative estimate.
C     DAMP for stabilization of Newton's method.
C     MXxITR is the maximum number of iterations.
C     DELDMX is the maximum step length.
C
      DATA       EPS    / 1.E-1 /  DAMP /1.E-4 /  MXxITR / 100 /
     .           DELDMX / 5.E-2 /
      EXTERNAL   FNCT
C>>>
      DATA NZERO /0/
C>>>
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      X     = XY(1)
      Y     = XY(2)
      DELX  = DXY(1)
      DELY  = DXY(2)
      XS    = X
      YS    = Y
      DELXS = DELX
      DELYS = DELY
      XSKL  = SKL(1)
      YSKL  = SKL(2)
      DTXMN = 1.E-4*ABS(XSKL)
      DTYMN = 1.E-4*ABS(YSKL)
      DTXMX = 1.E-2*ABS(XSKL)
      DTYMX = 1.E-2*ABS(YSKL)
      ADLF  = 1.E20
C--
C  Follow gradient vector to the point where FNCT(X,Y) = FCT.
C--
      NTRY  = 0
      DMP   = DAMP
1     NTRY  = NTRY + 1
      F     = FNCT(X,Y)
      DELF  = FCT - F
      ADLFS = ADLF
      ADLF  = ABS(DELF)
      DTX   = MAX(EPS*ABS(DELX),DTXMN)
      DTY   = MAX(EPS*ABS(DELY),DTYMN)
      DTX   = MIN(DTX,DTXMX)
      DTY   = MIN(DTY,DTYMX)
      DFX   = (FNCT(X+DTX,Y)-F)*XSKL/DTX
      DFY   = (FNCT(X,Y+DTY)-F)*YSKL/DTY
      IF ( ADLF.LE.FTEST )                       GO TO 10
C--
C  Decide next step.
C--
      ADF   = SQRT(DFX**2+DFY**2)
C>>> 91/05/30
      IF ( ADF.LE.1.D-20 ) THEN
         IF ( NZERO.LT.50 ) THEN
            NZERO = NZERO + 1
            PRINT *, ' >>>> ADF ZERO  ADF = ', ADF
         ENDIF
         X    = XS
         Y    = YS
         DELX = DELXS/2
         DELY = DELYS/2
         IF ( NTRY.LT.MXxITR )                   GO TO 1
                                                 GO TO 9
      ENDIF
C>>> 91/06/03
      IF ( ADLF.LT.ADLFS ) THEN
         DMP   = DMP/10
         XS    = X
         YS    = Y
         DELXS = DELX
         DELYS = DELY
C>>>
      ELSE
         DMP = DMP*10
      ENDIF
C>>>
      DELD  = SIGN(MIN(ADLF/ADF/(1+DMP),DELDMX),DELF)
      DELX  = DELD*XSKL*DFX/ADF
      DELY  = DELD*YSKL*DFY/ADF
C--
      X     = X + DELX
      Y     = Y + DELY
C--
      IF ( NTRY.LT.MXxITR )                      GO TO 1
9     PRINT *, ' Failed to get crossing point.'
CX    PRINT *, ' NTRY = ', NTRY, ' reached MXxITR in UGRADX.'
CX    PRINT *, '     ADLF = ', ADLF, ' FTEST = ', FTEST
CX    PRINT *, ' Restore XY, DXY.'
      XY(1)  = XS
      XY(2)  = YS
      DXY(1) = DELXS
      DXY(2) = DELYS
      NTRY   = - NTRY
C--
C  Now save (X,Y), (DELX,DELY), and (DFX,DFY)
C--
10    XY(1)  = X
      XY(2)  = Y
      DXY(1) = DELX
      DXY(2) = DELY
      DF(1)  = DFX/XSKL
      DF(2)  = DFY/YSKL
C--
C  That's it.
C--
      RETURN
      END
