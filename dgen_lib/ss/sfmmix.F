C----------------------------------------------------------------------C
      SUBROUTINE SFMMIX( WM, SIN2W,
     &                   HM1, TANB, ATERM,
     &                   SFTMSL, SFTMSR, FRMMAS, TYPFRM,
     &                   SFMAS, PHIF, THTF,
     &                   IRET  )
C----------------------------------------------------------------------C
C=======================================================================
C   Diagonalize mass matrices of sfermions
C   and then calculate some constants such as couplings etc.
C  Ref.   Hikasa-san's note
C
C
C <<< Inputs >>>  R*4 if not specified
C ( standard model constants )
C     1) WM   :    W-boson mass
C     2) SIN2W: weak mixing angle
C ( free parameters for SUSY )
C     1) HM1    :  the higgsino mass term
C     2) TANB   : the ratio of the v.e.v. of two higgs doublets
C       (  =V2/V1 )
C     3) ATERM(complex*8)  :  trilinear term A
C     4) SFTMSL,SFTMSR : soft scaler mass term of Left/Right-handed
C     5) FRMMAS : fermion mass
C     6) TYPFRM(I*4) : fermion type ( 0:neutrino, 1:charged lepton,
C                                  2:up-type quark, 3:down-type quark )
C
C    IPRNT ; print out flag for debug etc.
C
C <<< outputs >>>
C REAL*8  SFMAS(i) mass eigen values (i=1:lighter)
C REAL*8  PHIF     phase to make mass matrix real
C REAL*8  THTF     mixing angle
C     IRET < 0 : error
C     IRET > 0 : OK
C=============
C=======================================================================
 
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER ( PI = 3.1415927D0 )
      INTEGER*4 TYPFRM
      REAL*8    SFMAS(2), SFMAS2(2)
      COMPLEX*16 ATERM, AM2LR
 
C-----------------------------------------------------------------------
C check input (TYPFRM) & set T3L,charge
      IF( TYPFRM.EQ.0 ) THEN
C        neutrino
        T3L = 0.5D0
        QCHRG = 0
      ELSEIF( TYPFRM.EQ.1 ) THEN
C        charged lepton
        T3L = -0.5D0
        QCHRG = -1
      ELSEIF( TYPFRM.EQ.2 ) THEN
C        up-type quark
        T3L = 0.5D0
        QCHRG = 2/3.D0
      ELSEIF( TYPFRM.EQ.3 ) THEN
C        down-type quark
        T3L = -0.5D0
        QCHRG = -1/3.D0
      ELSE
C       WRITE(6,*) ' Error(SFMMIX):TYPFRM=',TYPFRM
        IRET = -1
        RETURN
      ENDIF
 
      COS2W = 1-SIN2W
      TAN2W = SIN2W/COS2W
      SINW  = SQRT(SIN2W)
      COSW  = SQRT(COS2W)
      TANW  = SQRT(TAN2W)
      ZM    = WM/COSW
      WM2   = WM**2
      ZM2   = ZM**2
 
C  0 < beta < pi
      COSB  = SQRT( 1./(TANB**2+1.) )
      IF( TANB .LT. 0 ) COSB = -COSB
      SINB  = SQRT((1-COSB)*(1+COSB))
      COSDB = (COSB-SINB)*(COSB+SINB)
      SINDB = 2*SINB*COSB
      SIN2B = SINB**2
      COS2B = COSB**2
 
C off-diagonal element
      IF( MOD(TYPFRM,2).EQ.0 ) THEN
        AM2LR = -FRMMAS*( CONJG(ATERM) +HM1/TANB )
      ELSE
        AM2LR = -FRMMAS*( CONJG(ATERM) +HM1*TANB )
      ENDIF
 
      SFTML2 = SFTMSL**2
      SFTMR2 = SFTMSR**2
      FRMMS2 = FRMMAS**2
 
      SM2LR = ABS(AM2LR)**2
      ASUM1 = SFTML2 + FRMMS2
      SM2FL = ASUM1 + ZM2*COSDB*( T3L-QCHRG*SIN2W)
      SM2FR = ASUM1 + ZM2*COSDB*QCHRG*SIN2W
 
      CM2LR = REAL(AM2LR)
      SM2LR = IMAG(AM2LR)
      PHIF  = ATAN2( SM2LR, CM2LR )
 
      ASUM1  = SFTML2 +SFTMR2
      ASQRT  = SQRT( (SFTML2-SFTMR2)**2 + 4*SM2LR**2 )
      SFMAS2(1) = ( ASUM1 - ASQRT )/2
      SFMAS2(2) = ( ASUM1 + ASQRT )/2
C>>> 96/08/02  Bug fixup due to nojirim
C      DO 10 I=1,2
C10    SFMAS(I) = SQRT(SFMAS(I))
      DO 10 I=1,2
10    SFMAS(I) = SQRT(ABS(SFMAS2(I)))
      IF(SFMAS2(1).LT.0) THEN 
        IRET=-1
        RETURN
      ENDIF
C>>>
 
      TANTHF = SM2LR/(SFMAS2(1) - SFTMR2 )
      THTF   = ATAN2( TANTHF, 1.D0 )
 
      IRET = 1
 
      RETURN
      END
