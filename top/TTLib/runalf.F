***********************************************************************
*  REAL*8 FUNCTION ALFMU( MU, ALFZ )
*
*	This function subroutine calculates the running coupling 
*	constant alpha_s^MSbar(mu).
*	Input parameter is alpha_s^MSbar(mZ).
***********************************************************************

      real*8  function  alfmu( mu, alfz )


      implicit none

      real*8  alfz, mu

      real*8  CF, CA, TF, mZ, nf

      real*8  runalf

      real*8  pi, b0, b1, Lam_msbar
      common/run/pi, b0, b1, Lam_msbar


      pi = 3.14159265
      CF = 4./3.
      CA = 3.
      TF = 1./2.
      nf = 5.
      mZ = 91.2


      b0 = 11./6.*CA -  1./3.*nf
      b1 = 17./12.*CA**2 - 5./6.*CA*nf*TF - 1./2.*CF*nf*TF

      Lam_msbar = mZ*exp( - pi/b0/alfz
     -		  + b1/b0**2*log( 2./b0*(pi/alfz+b1/b0) ) )


      alfmu = runalf(mu)

C    OUTPUT

      write(*,*) 'Lamda_MSbar =', Lam_msbar
C      write(*,*) 'alpha_s(',real(mu),') =', alfmu

      return
      end


*****************************************************************

      real*8  function runalf(mu)


      implicit none

      real*8  mu

      real*8  pi, b0, b1, Lam_msbar
      common/run/pi, b0, b1, Lam_msbar

      real*8  g, f, fd, x, x0, x1


      g(x) = mu*exp( - pi/b0/x 
     -	             + b1/b0**2*log( 2./b0*(pi/x+b1/b0) ) )

      f(x) = g(x) - Lam_msbar

      fd(x) = g(x)*( pi/b0/x**2 - b1/b0**2
     -		     * ( 2.*pi/b0/x**2 )/( 2./b0*(pi/x+b1/b0) ) )


      x0 = (pi/b0)
     -   /( log(mu/Lam_msbar) + (b1/b0**2)*log(log(mu/Lam_msbar)) )

      do while ( abs(f(x0)) .gt. 1.d-4 )

	  x1 = x0 - f(x0)/fd(x0)
	  x0 = x1

      end do

      runalf = x0

      return
      end
