**********************************************************************
*
* Subroutine QCDMOD_P(ALFS,AMT,GMT,AMH,BTH,RS,Q12,Q22,PQCD)
*
* (Purpose)
*    Calculates QCD modification factor to the ttbar axial-vector
*    vertex determined from the EW amplitude for e+e- --> bW+ bbarW-.
* (Inputs)
*    ALFS   : (R*4) : alpha_s.
*    AMT    : (R*4) : top mass.
*    GMT    : (R*4) : top width.
*    AMH    : (R*4) : Higgs mass.
*    BTH    : (R*4) : top's Yukawa coupling normalized by SM.
*    RS     : (R*4) : sqrt(s).
*    Q12    : (R*4) : virtual top mass squared
*    Q22    : (R*4) : virtual anti-top mass squared.
* (Output)
*    PQCD   : (C*8): (Q1,Q2)-dependent QCD modification to the ttbar 
*                  : axial-vector coupling calculated without QCD.
* (Update Record)
*   99/7/26  Y.Sumino        Original version.
* (Usage)
*   Call "tabulate_GreenP" before using this subroutine.
*
**********************************************************************
 
      SUBROUTINE QCDMOD_P(ALFS,AMT,GMT,AMH,BTH,RS,Q12,Q22,PQCD)
 
      implicit none

      REAL*4      ALFS, AMT, GMT, AMH, BTH, RS, Q12, Q22
      COMPLEX*8   PQCD

      real*4  E, pt

      complex*8  G_Pwave

      real*4  pi
      complex*8  imag
      integer*4 ncall
      data ncall   /0/
C--
C  Initialization.
C--
      if ( ncall.eq.0 ) then
         ncall = 1
         print *, 'QCDMOD_P now tabulates GQ_P'
         call tabulate_greenp(alfs,amt,gmt)
         pi = acos(-1.d0)
         imag = ( 0.d0, 1.d0 )
      endif

      E = RS - 2.*AMT

      pt = sqrt( ( (RS**2-Q22+Q12)/(2.*RS) )**2 - Q12 )

      if ( ( E .gt. -6. ) .and. ( E .lt. 10. ) 
     -     .and. ( pt .lt. 100. ) ) then

         PQCD = - ( 1. - 4.*ALFS/3./pi )
     -        *( E + imag*GMT - pt**2/AMT )*G_Pwave(E,pt)

      else

         PQCD = 1.

      endif

      return

      end

***********************************************************************

      COMPLEX*8  FUNCTION  G_PWAVE(ENER,P)

      implicit none


      real*4  ener, p

      integer  n_ener, n_p

      real*4  dx, dy, dp

      complex*8  array_GP(0:40,0:800)
      common/tabulateGP/array_GP

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

      integer     n_green
      common/size/n_green


      n_ener = int( 2.*(ener+6.) )
      dx = ( ener + 6. - real(n_ener)/2. )/0.5

      dp = 2.*Lambda/n_green
      n_p = int( p/dp )
      dy = ( p - real(n_p)*dp )/dp


      G_Pwave = (1.-dx)*(1.-dy)*array_GP( n_ener,   n_p   )
     -          + dx*(1.-dy)   *array_GP( n_ener+1, n_p   )
     -          + (1.-dx)*dy   *array_GP( n_ener,   n_p+1 )
     -          + dx*dy        *array_GP( n_ener+1, n_p+1 )

C      G_Pwave = (
C     -    (3.-2.*dx-2.*dy)*array_GP( n_ener,   n_p   )
C     -  + (2.*dx-2.*dy+1.)*array_GP( n_ener+1, n_p   )
C     -  + (2.*dy-2.*dx+1.)*array_GP( n_ener,   n_p+1 )
C     -  + (2.*dx+2.*dy-1.)*array_GP( n_ener+1, n_p+1 )
C     -           )/4.

C      G_Pwave = (1.-dx-dy)*array_GP( n_ener,   n_p   )
C     -         +    dx    *array_GP( n_ener+1, n_p   )
C     -         +    dy    *array_GP( n_ener,   n_p+1 )

C      write(*,*)  'dx, dy', dx, dy
C      write(*,*) array_GP( n_ener,   n_p   )
C      write(*,*) array_GP( n_ener+1, n_p   )
C      write(*,*) array_GP( n_ener,   n_p+1 )
C      write(*,*) array_GP( n_ener+1, n_p+1 )

      return

      end


