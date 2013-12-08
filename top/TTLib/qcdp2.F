*********************************************************************
*
*  SUBROUTINE QCDP1( E, DR, EPS, STEP1 )
*
*  This subroutine solves the P-wave Schroedinger equation:
*
*        dy/dx = z
*        dz/dx = ( -eff_mt*(E+iGamma0/2-V(x)) + 2/x^2 )*y
*
*  using the Lunge-Kutta method.
*
*  It determines the two independent solutions
*  y1(x) and y2(x), which are determined by the conditions at
*  x=0.01 GeV^-1 :
*
*        y1 = ( 4.d-3, -7.d-4 ), z1 = ( 4.d-2, -2.d-2 )
*        y2 = ( 2. d8, -3. d7 ), z2 = ( 2. d9, -1. d9 ) .
*
*  y1(x) and y2(x) are given in common/PFRSTL/... for x>0.01,
*  and in common/PFRSTS/... for x<0.01 .
*
*  Input:
*       E     : energy = sqrt(s) - 2*m_t  ;  it will be stored in
*                                common/enerP/... to be used in the
*                                function subroutine  g_P .
*       dr    : integration step used in the region x>0.01
*       eps   : criterion for the convergence of B;
*               integration continued while  |B_old/B_new-1| > eps
*
*  Output: x_L(num), y1_L(num), y2_L(num), and num_L
*                               are stored in common/PFRSTL/...
*          x_S(num), y1_S(num), y2_S(num), and num_S
*                               are stored in common/PFRSTS/...
*
*       step1 : number of steps required in x >= 0.01 to obtain
*               convergent  B = - lim_{x->infty} y1/y2
*       x_L, x_S   : x < 0.01 and x > 0.01, respectively
*       y1_L, y1_S : solution y1 corresponding to x_L and x_S, rsptvly
*       y2_L, y2_S : solution y2 corresponding to x_L and x_S, rsptvly
*       num_L : number of points for x >= 0.01
*       num_S : number of points for x < 0.01
*
*********************************************************************


      subroutine QCD_P1( E, dr, eps, step1 )


      implicit none


C    INPUT VARIABLES

      real*8  E, dr, eps

C    OUTPUT VARIABLES

#if 0
      real*8  x_L(5000), x_S(5000)
      complex*16  y1_L(5000), y2_L(5000)
      complex*16  y1_S(5000), y2_S(5000)
      integer num_L, num_S
      common/PFIRSTL/ x_L, y1_L, y2_L, num_L
      common/PFIRSTS/ x_S, y1_S, y2_S, num_S
#else
#include "pfirst.inc"
#endif

      integer  step1


C    OTHER VARIABLES

      real*8  enerP
      common/enerP/enerP

      real*8  x, dx, criterion
      integer maxstps, n, clock, k, num

      complex*16 y1, y2, z1, z2, 
     -        dy1, dz1, dy2, dz2,
     -        B, B_new, r1, r2, r1_new, r2_new


C    MAXIMUM NUMBER OF STEPS ALLOWED
      parameter( maxstps = 400000 )
C      parameter( maxstps = 40000 )
#if 1
      integer    n_mx
      parameter( n_mx = 45000 )
#endif

C      write(*,*)
      write(*,*) 'P-wave: E=', E


      enerP = E


C  FOR r>=0.01 *****************************

C    CONDITIONS AT X = 0.01

      x = 0.01

      y1 = (4.d-3,-7.d-4)
      z1 = (4.d-2,-2.d-2)

      y2 = (2.d8,-3.d7)
      z2 = (2.d9,-1.d9)

      x_L(1) = x
      y1_L(1) = y1
      y2_L(1) = y2


C    DO LOOP FOR X > 0.01  CONTINUED WHILE  ( CRITERION > EPS )

      n = 0
      B = 0.
      criterion = 100.
      dx = dr
      num = 1

#if 0
      do while ( ( criterion .gt. eps ) .or.
     -                       ( n .lt. 500 ) )
#else
      do while ( ( criterion .gt. eps ) .or.
     -                       ( n .lt. n_mx ) )
#endif

        do k = 1, 2

          do clock = 1, 20

              call  lungekutta_P(x,y1,z1,dx,dy1,dz1)
              call  lungekutta_P(x,y2,z2,dx,dy2,dz2)

              y1 = y1 + dy1
              z1 = z1 + dz1

              y2 = y2 + dy2
              z2 = z2 + dz2

	      x = x + dx

              n = n + 1
#if 1
              if (n.ge.n_mx) goto 90
#endif
          end do

  	  num = num + 1

          x_L(num) = x
          y1_L(num) = y1
          y2_L(num) = y2

C          write(*,*) x_L(num)
C          write(*,*) y1_L(num)
C          write(*,*) y2_L(num)

        end do
#if 1
90      continue
#endif
C    CALCULATE CRITERION

        B_new = - y1/y2

        criterion = abs(B/B_new-1.)

        B = B_new

        if ( n .gt. maxstps ) then
            write(*,*)  'integration steps exceeded',maxstps
	    write(*,*)  criterion
            goto 100
        end if

      end do

C    END OF DO LOOP

      num_L = num
      step1 = num

C      write(*,*) 'r>=0.01 completed: step1, B, B_new, criterion'
C      write(*,*)  step1
C      write(*,*)  B, B_new
C      write(*,*)  criterion


C  FOR r<0.01 *****************************

C    CONDITIONS AT X = 0.01

      x = 0.01

      y1 = (4.d-3,-7.d-4)
      z1 = (4.d-2,-2.d-2)

      y2 = (2.d8,-3.d7)
      z2 = (2.d9,-1.d9)


C    DO LOOP FOR X < 0.01  CONTINUED WHILE  ( CRITERION > 1.d-3 )

      n = 0
      r1 = 0.
      r2 = 0.
      criterion = 100.
      num = 0

      do while ( criterion .gt. 1.d-3 )

        dx = -x*.02

        do k = 1, 2

          do clock = 1, 10

              call  lungekutta_P(x,y1,z1,dx,dy1,dz1)
              call  lungekutta_P(x,y2,z2,dx,dy2,dz2)

              y1 = y1 + dy1
              z1 = z1 + dz1

              y2 = y2 + dy2
              z2 = z2 + dz2

	      x = x + dx

              n = n + 1

          end do

  	  num = num + 1

          x_S(num) = x
          y1_S(num) = y1
          y2_S(num) = y2

        end do

        r1_new = x * y1
        r2_new = x * y2

        criterion = abs(r1/r1_new-1.)+abs(r2/r2_new-1.)

        r1 = r1_new
        r2 = r2_new

        if ( n .gt. maxstps ) then
            write(*,*)  'integration steps exceeded',maxstps
	    write(*,*)  criterion
            goto 100
        end if

      end do


C    END OF DO LOOP

      num_S = num

C      write(*,*) 'r<0.01 completed: n, r1, r2, criterion'
C      write(*,*) n, r1, r2
C      write(*,*) criterion
C      write(*,*) 'x =', x
C      write(*,*) 'steps =', num


 100  continue


      return

      end


*******************************************************************
*  SUBROUTINE LUNGEKUTTA_P(X,Y,Z,DX,DY,DZ)
*
*  This subroutine solves the differential equation
*
*     dy/dx = f_P(x,y,z),  dz/dx = g_P(x,y,z)
*
*  using the Lunge-Kutta method.  x is the real variable,
*  while y, z, f_P, g_P are complex.
*
*  Inputs:  x, y, z, dx
*  Outputs: dy, dz
*
********************************************************

      subroutine lungekutta_P(x,y,z,dx,dy,dz)


      implicit none


      real*8 x, dx, half_dx

      complex*16 y, z, dy, dz, d1, d2, d3, d4, 
     -        e1, e2, e3, e4, f_P, g_P



      half_dx = .5*dx


      d1 = dx*f_P( x, y, z )
      e1 = dx*g_P( x, y, z )

      d2 = dx*f_P( x+half_dx, y+.5*d1, z+.5*e1 )
      e2 = dx*g_P( x+half_dx, y+.5*d1, z+.5*e1 )

      d3 = dx*f_P( x+half_dx, y+.5*d2, z+.5*e2 )
      e3 = dx*g_P( x+half_dx, y+.5*d2, z+.5*e2 )

      d4 = dx*f_P( x+dx, y+d3, z+e3 )
      e4 = dx*g_P( x+dx, y+d3, z+e3 )

      dy = ( d1 + 2.*d2 + 2.*d3 + d4 )/6.
      dz = ( e1 + 2.*e2 + 2.*e3 + e4 )/6.

      return
      end


*********************************************************
*  COMPLEX*16 FUNCTION  F_P(X,Y,Z)
*********************************************************

      complex*16 function f_P(x,y,z)

      implicit none


      real*8 x

      complex*16 y, z


      f_P = z


      end


*********************************************************************
*  COMPLEX*16 FUNCTION  G_P(X,Y,Z)
*
*  Inputs: 
*	x, y, z,
*	enerP : E = sqrt(s) - 2m_t ; stored in  common/enerP/...
*	mt, eff_mt, imag;  stored in common/param/...
*  Output:
*	g_P = ( - eff_mt*( E + i*Gamma0/2 - V(x) ) + 2/x**2 )*y
*
**********************************************************************

      complex*16 function g_P(x,y,z)

      implicit none


C    INPUT VARIABLES

      real*8  enerP
      common/enerP/enerP


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



C    OTHER VARIABLES

C      real*8 x, Vc_tab_118_3
      real*8 x, V
      complex*16 y, z


C    OUTPUT

C      g_P = ( - eff_mt*( enerP + imag*Gammat - Vc_tab_118_3(x) )
      g_P = ( - eff_mt*( enerP + imag*Gammat - V(x) )
     -	                                               + 2./x**2)*y


      end


*****************************************************************
*
*  SUBROUTINE QCDP2( STEP2 )
*
*  This subroutine calculates the factor B_P from the
*  asymptotic behavior of y1 and y2.
*  Then determine f(x) which is the solution to the
*  homogeneous P-wave Schroedinger equation with the
*  boundary condition
*
*       f(0) = 1,    f(x) -> 0  as x->infty .
*
*  This solution is given by  f(x) = A [ y1(x) + B_P*y2(x) ].
*
*  Inputs:
*       x_L, y1_L, y2_L : two independent solutions at x>=0.01
*                         ; stored in  common/PFIRSTL/...
*       x_S, y1_S, y2_S : two independent solutions at x<0.01
*                         ; stored in  common/PFIRSTS/...
*       num_L, num_S    : number of points for the solutions,
*                         for x>=0.01 and x<0.01, respctvly.
*  Outputs:
*       dst_P(num), f(num)  : solution satisfying the above boundary
*                         conditions ; stored in  common/PSECOND/...
*       num_max         : number of points for f(x) ;
*                               also stored in common/PSECOND/...
*       step2 : number of steps required for sufficiently convergent
*               f(x) for x>=0.01
*
*********************************************************************

      subroutine  QCD_P2( step2 )

      implicit none


C    INPUT VARIABLES

#if 0
      real*8  x_L(5000), x_S(5000)
      complex*16  y1_L(5000), y2_L(5000)
      complex*16  y1_S(5000), y2_S(5000)
      integer num_L, num_S
      common/PFIRSTL/ x_L, y1_L, y2_L, num_L
      common/PFIRSTS/ x_S, y1_S, y2_S, num_S
#else
#include "pfirst.inc"
#endif


C    OUTPUT VARIABLES

#if 0
      real  dst_P(5000)
      complex  f(5000)
      integer num_max
      common/PSECOND/ dst_P, f, num_max
#else
#include "psecond.inc"
#endif

      integer  step2


C    OTHER VARIABLES

      real*8 norm, norm_max, criterion

      complex*16 A, B_P, B_tmp, T

      integer num, k

      integer n


C    SET CRITERION FOR TRANCATING THE GREEN'S FUNCTION F(X)

      parameter ( criterion = 1.d-5 )



C    DETERMINATION OF B_P = - LIM_{X->INFTY} Y1(X)/Y2(X)

      T = 0.
      do num = num_L-4, num_L

          B_tmp = - y1_L(num)/y2_L(num)
          T = T + B_tmp

      end do

      B_P = T/5.


C    DETERMINATION OF A = LIM_{X->0} 1/{ X*[ Y1 + B_P*Y2 ] }

      A = 1./x_S(num_S)/( y1_S(num_S) + B_P*y2_S(num_S) )


C    DETERMINATION OF F(X) FOR X<0.01

      num = 0
      do k = num_S, 1, -1

          num = num + 1
 	  dst_P(num) = x_S(k)
          f(num) = A*( y1_S(k) + B_P*y2_S(k) )

C          write(50,*) dst_P(num), abs(f(num))**2

      end do



C    DETERMINATION OF F(X) FOR X>=0.01

      k = 0
      norm_max = 0.
      do while ( k .lt. num_L )

	  k = k + 1
	  num = num + 1
	  dst_P(num) = x_L(k)
          f(num) = A*( y1_L(k) + B_P*y2_L(k) )

          norm = dst_P(num)**2 * abs(f(num))**2
          norm_max = max( norm, norm_max )

C          write(50,*) dst_P(num), norm

          if ( norm .lt. criterion*norm_max ) goto 200

      end do

      write(*,*) 'f(r) does not become sufficiently small.'


 200  continue

      step2 = k
      num_max = int((num-1)/2)*2 + 1

C      write(*,*) 'total number of points, used number of points'
C      write(*,*) num_S+num_L, num_max


      return

      end

      
*********************************************************************
*
*  SUBROUTINE QCDP3
*
*  This subroutine calculates the Fourier Transform
*  of F^l(x;E) :
*
*                                  _infty
*                         eff_mt  /
*         tilde{F}(p;E) = ------ /  dr  f(r)*delfn(pr),
*                          p^2 _/ 0
*
*  where the Fourier kernel delfn(x) is defined as
*
*         delfn(x) = (sin(x)-x*cos(x))/x .
*
*  tilde{F}(p;E) is given at ( n_green + 1 ) points for 0 =< p =< 2*Lambda
*
*  Inputs:
*       dst_P(num), f(num) : function f(r) ;
*       num_max : number of points of f(r) ;
*                               stored in common/PSECOND/...
*       eff_mt  : top quark mass ;
*       Lambda  : maximum momentum ;
*                               stored in common/param/..
*  Outputs:
*       p_P(k), tilf(k)  : Green's function in the momentum space
*                       = tilde{F}(p;E)  ;
*                               stored in common/PTHIRD/...
*
*  Remark:  We set p_P(0)=0, tilf(0)=0 .
*
*********************************************************************

      subroutine QCD_P3( n_E )


      implicit none


C    INPUT VARIABLES

#if 0
      real  dst_P(5000)
      complex  f(5000)
      integer num_max
      common/PSECOND/ dst_P, f, num_max
#else
#include "psecond.inc"
#endif

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

C    OUTPUT VARIABLES

      real  p_P(0:5000)
      complex  tilf(0:5000)
      common/PTHIRD/ p_P, tilf

#if 0
      complex*8  array_GP(0:40,0:800)
      common/tabulateGP/array_GP
#else
#include "tabgp.inc"
#endif

C    OTHER VARIABLES

      integer n, k, n_E
      real q, r1, r2, r3, criterion
      complex f1, f2, f3, integrand

      real delfn


      parameter( criterion = 1.e-3 )


C    FOR p = 0

      p_P(0) = 0.
      tilf(0) = 0.


C    DO LOOP FOR  0 < p =< 2*Lambda in n_green steps


      do n=1, n_green

          q = 2.*Lambda*real(n)/n_green

          integrand = 0.

	  r1 = dst_P(1)
          f1 = f(1)

          do k = 2, num_max-1, 2

	      r2 = dst_P(k)
	      f2 = f(k)

	      r3 = dst_P(k+1)
	      f3 = f(k+1)

	      if ( ((r3-r2)/(r2-r1)-1.) .gt. criterion )
     - 		write(*,*)  'integration steps not equal'

              integrand = integrand 
     -     + ( f1*delfn(q*r1) + 4.*f2*delfn(q*r2) 
     -                  + f3*delfn(q*r3) )* (r3-r1)

              r1 = r3
              f1 = f3

          end do

	  p_P(n) = q
          tilf(n) = eff_mt/q**2 /6.*integrand

          array_GP( n_E, n ) = tilf(n)

C          write(*,*) 'n_E, n:', n_E, n
C          write (40,*)  q, q**2*abs(tilf(n))**2 


      end do

C    DO LOOP COMPLETED

      return

      end

*******************************************************
*   real function delfn(x)
*******************************************************

      real function delfn(x)

      implicit none


      real x

      delfn = ( sin(x) - x*cos(x) )/x

      end

************************************************************





