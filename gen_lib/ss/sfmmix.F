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
C ( free parameters for SUSY )(
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
C REAL*4  SFMAS(i) mass eigen values (i=1:lighter)
C REAL*4  PHIF     phase to make mass matrix real
C REAL*4  THTF     mixing angle
C     IRET < 0 : error
C     IRET > 0 : OK
C=============
C=======================================================================
 
      PARAMETER ( PI = 3.1415927 )
      INTEGER*4 TYPFRM
      REAL*4    SFMAS(2), SFMAS2(2)
      REAL*4 ATERM, AM2LR
 
C-----------------------------------------------------------------------
C check input (TYPFRM) & set T3L,charge
      IF( TYPFRM.EQ.0 ) THEN
C        neutrino
        T3L = 0.5
        QCHRG = 0.
      ELSEIF( TYPFRM.EQ.1 ) THEN
C        charged lepton
        T3L = -0.5
        QCHRG = -1.
      ELSEIF( TYPFRM.EQ.2 ) THEN
C        up-type quark
        T3L = 0.5
        QCHRG = 2./3.
      ELSEIF( TYPFRM.EQ.3 ) THEN
C        down-type quark
        T3L = -0.5
        QCHRG = -1./3.
      ELSE
C       WRITE(6,*) ' Error(SFMMIX):TYPFRM=',TYPFRM
        IRET = -1
        RETURN
      ENDIF
 
      COS2W = 1.-SIN2W
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
      SINB  = SQRT((1.-COSB)*(1.+COSB))
      COSDB = (COSB-SINB)*(COSB+SINB)
      SINDB = 2.*SINB*COSB
      SIN2B = SINB**2
      COS2B = COSB**2
 
C off-diagonal element
      IF( MOD(TYPFRM,2).EQ.0 ) THEN
        AM2LR = -FRMMAS*(ATERM +HM1/TANB )
      ELSE
        AM2LR = -FRMMAS*( ATERM +HM1*TANB )
      ENDIF
 
      SFTML2 = SFTMSL**2
      SFTMR2 = SFTMSR**2
      FRMMS2 = FRMMAS**2
 
      SM2LR = AM2LR
C      SM2LR=0.
      ASUM1 = SFTML2 + FRMMS2
      SM2FL = ASUM1 + ZM2*COSDB*( T3L-QCHRG*SIN2W)
      ASUM2 = SFTMR2 + FRMMS2
      SM2FR = ASUM2 + ZM2*COSDB*QCHRG*SIN2W
      DL=ZM2*COSDB*( T3L-QCHRG*SIN2W)
      DR=ZM2*COSDB*QCHRG*SIN2W
C      WRITE(*,*) 'DTERMS', DL,DR, 'COSDB=', COSDB
      ASUM1  = SM2FR +SM2FL
      ASQRT  = SQRT( (SM2FR-SM2FL )**2 + 4*SM2LR**2 )
      SFMAS2(1) = ( ASUM1 - ASQRT )/2.
      SFMAS2(2) = ( ASUM1 + ASQRT )/2.
      DO 10 I=1,2
10    SFMAS(I) = SQRT(ABS(SFMAS2(I)))
C      write(*,*) 'from sfmmix'
C      write(*,*) sm2lr, sfmas2(1), sm2fr
C
C
      if(sm2lr.gt.0.0001*abs(sm2fl-sm2fr)) then 
      TANTHF = SM2LR/(SFMAS2(1) - Sm2FR )
      THTF   = ATAN2( TANTHF, 1. )
      elseif(sm2fr.gt.sm2fl)  then 
      TANTHF = SM2LR/(SFMAS2(1) - Sm2FR )
      THTF   = ATAN2( TANTHF, 1. )
      elseif(sm2fr.le.sm2fl) then 
      thtf  = atan2((sm2fl-sm2fr), sm2lr)
      endif
C
      IF(SFMAS2(1).LT.0) THEN 
        IRET=-1
        RETURN
      ENDIF
      IRET = 1
      RETURN
      END





