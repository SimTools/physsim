***********************************************************************
*     TABULATE_GREENP.F
*
*  This subroutine tabulates the P-wave Green function.
*  
* (Inputs)
*    ALFS   : (R*4) : alpha_s.
*    AMT    : (R*4) : top mass.
*    GMT    : (R*4) : top width.
* (Outputs)
*    ARRAY_GP(-12:20,0:500) : (R*4) : table of P-wave Green function
*     stored in   common/tabulateGP/
***********************************************************************

      SUBROUTINE  TABULATE_GREENP(ALFS,AMT,GMT)

      implicit none


C   INPUT VARIABLES

      real*4  ALFS, AMT, GMT

C   OTHER VARIABLES

      real*8  mt, alfz, E, y_cut,
     -        drS, drP, epsS, epsP
      real*8  r_wt, Gammat, Lambda, CF, mb, mZ, 
     -        tileta1, tileta2, kappa
      real*8  pi
      real*8  alfp
      complex*16  eff_mt
      complex*16  imag

      common/param/mt, alfz, E, y_cut, 
     -             drS, drP, epsS, epsP,
     -             r_wt, Gammat, Lambda, CF, mb, mZ, 
     -             eff_mt, tileta1, tileta2, kappa,
     -		   pi, imag 

      real*8  alft, r_QCD

      common/QCDpar/alft, r_QCD


C    SET UP PARAMETERS
      call setparam(ALFS,AMT,GMT)

C    TABULATE QCD POTENTIAL
      call setvpr(alfz,mt,mb,mZ,Gammat,alft,alfp,r_QCD)
c      call PREPARE_VcRG

C    MAIN ROUTINE
      call routine


      end

*******************************************************************
*******************************************************************
*  SUBROUTINE SETPARAM
*
*  This subroutine prepares the parameter data set necessary
*  for the calculation of cross sections.
*
*  Inputs:
*    ALFS   : (R*4) : alpha_s.
*    AMT    : (R*4) : top mass.
*    GMT    : (R*4) : top width.
*  Output:    stored in  common/param/...
*	   mt, alfz, E, y_cut, drS, drP,
*	   r_wt   : mW**2/mt**2
*	   Gammat : top quark on-shell width
*          Lambda : cut-off of the top quark three-momentum
*          CF     : color factor = 4/3
*	   eff_mt : complex effective mass related to the running 
*		    toponium width [ = mt/( 1 + i*kap ) ]
*	   tileta1, tileta2 : coefficients of the running toponium 
*		    width including the order alpha_s correction
*	   kappa  : coupling related to the axial-vertex
*		    for t->bW  [ = (1-2r)/(1+2r) ]
*	   mb, mZ, pi, imag
*****************************************************************

      subroutine  setparam(ALFS,AMT,GMT)


      implicit none


C   INPUT VARIABLES

      real*4  ALFS, AMT, GMT

C   OTHER VARIABLES

      real*8  mt, alfz, E, y_cut,
     -        drS, drP, epsS, epsP
      real*8  r_wt, Gammat, Lambda, CF, mb, mZ, 
     -        tileta1, tileta2, kappa
      real*8  pi
      complex*16  eff_mt
      complex*16  imag

      common/param/mt, alfz, E, y_cut, 
     -             drS, drP, epsS, epsP,
     -             r_wt, Gammat, Lambda, CF, mb, mZ, 
     -             eff_mt, tileta1, tileta2, kappa,
     -		   pi, imag 

      real*8  alfsmt, alfsmuB, muB
      common/coupl/alfsmt, alfsmuB, muB

      real*8  alfmu

      real*8  Gammat0, GF, mW, radcorr, h, kap2

      real*8  width_withmb

      real*8  zero

      integer     n_green
      common/size/n_green


C    PARAMETER SETUP

      pi   = 3.1415926
      imag = (0.,1.)

      mb   = 5.
      mZ   = 91.2
      mW   = 80.6

      CF   = 4./3.

      GF   = 1.166E-5

      muB = 15.

      zero = 0.

      n_green = 200

C    INPUTS

      mt = AMT
      alfz = ALFS
      Gammat = GMT

      write(*,*) mt,alfz,Gammat
C      write(*,*) 'm_t, alpha_s(mZ)?:'
C      read(*,*) mt, alfz

C    SETTING UP DRS, DRP, EPSS, EPSP

      drS = .01/mt
      drP = .01/mt

      epsS = 1.d-10
      epsP = 1.d-10

C    CALCULATION OF COUPLING CONSTANT AT VARIOUS SCALES

      alfsmt  = alfmu(  mt, alfz )
      alfsmuB = alfmu( muB, alfz )

      write(*,*) 'alfsmt, alfsmuB'
      write(*,*) alfsmt, alfsmuB

C    CALCULATION OF LAMBDA = CUT-OFF

      Lambda = sqrt( (9*mt**2-mW**2)*(mt**2-mW**2) )/4./mt

C    CALCULATION OF RELICS

      r_wt = mW**2/mt**2
      kappa = ( 1. - 2.*r_wt )/( 1. + 2.*r_wt )

      tileta1 = 0.
      tileta2 = 0.

      kap2 = Gammat*tileta2/mt

      eff_mt = mt/( 1. + imag*kap2 )

      return

      end



************************************************************************
C      include 'runalf.for'

************************************************************************
*  SUBROUTINE PREPARE_VcRG
C      include 'qcdpot.for'
C      include 'Vc_tab_118_3.for'

*****************************************************************************
*  SUBROUTINE ROUTINE
*
*  This subroutine calculates the Green functions.
*  
*  Input: E	stored in common/param/...
*
*  Output: File = 'GreenS.dat'
*	   File = 'GreenP.dat'
*
*****************************************************************

      subroutine routine

      implicit none


      real*8  mt, alfz, E, y_cut,
     -        drS, drP, epsS, epsP
      real*8  r_wt, Gammat, Lambda, CF, mb, mZ, 
     -        tileta1, tileta2, kappa
      real*8  pi
      complex*16  eff_mt
      complex*16  imag

      common/param/mt, alfz, E, y_cut, 
     -             drS, drP, epsS, epsP,
     -             r_wt, Gammat, Lambda, CF, mb, mZ, 
     -             eff_mt, tileta1, tileta2, kappa,
     -		   pi, imag 

      integer  stepS1, stepS2, stepP1, stepP2,
     -	   stepS2_old, stepP2_old

      integer  n, n_E


C    START

      do n = -12, 20

          E = real(n)/2.
          n_E = n + 12

          call  QCD_P1(E,drP,epsP,stepP1)

          call  QCD_P2(stepP2)

          call  QCD_P3(n_E)

      enddo

      return
      end


****************************************************************
*  SUBROUTINE QCD_P1( E, DR, EPS, STEP1 )
*****************************************************************
*  SUBROUTINE QCD_P2( STEP2 )
**********************************************************************
*  SUBROUTINE QCD_P3

C      include 'qcdp2.for'


