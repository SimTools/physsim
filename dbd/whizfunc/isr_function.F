      subroutine isr_function(factor, x, eps, LLA_order)
c
c(function)
c  generate ISR spectrum a la Jadaha.
c  This code is extracted from Whizard/whizard-src/beams.f90
c(arguments)
c  factor [R*8,inout]
c  x      [R*8,inout]
c  eps    [R*8,in]
c  LLA_order [I*4, in]
c
c(Author)
c  Akiya Miyamoto, 2-June-2011, original version
c
c ========================================================
c      implicit none
      double precision factor, x, eps
      integer LLA_order

      double precision f
      double precision xx, xxbar, xbar, log_x, log_xbar, x_2

      double precision xmin, zeta3
      double precision pi, g1, g2, g3
      double precision Li2_ISR
      parameter (xmin=0.00714053329734592839549810275603D0)
      parameter (zeta3 = 1.20205690315959428539973816151D0)
      parameter (pi=3.14159265358979D0)
      parameter (g1=3.0D0/4.0D0)
      parameter (g2=(27.0D0-8.0D0*pi**2)/96.0D0) 
      parameter (g3=(27.0D0-24.0D0*pi**2+128.0D0*zeta3)/384.0D0)

c ================< entry point > =======================
      xx = x
      xxbar = 1 - xx
      if (xxbar < tiny(1.0D0)**eps) then
        xbar = 0
      else
        xbar = xxbar**(1.0D0/eps)
      end if
      x = 1 - xbar

c       print *,' isr_function ... xx=',xx,' x=',x,
c     >          ' eps=',eps,' LLA_Order=',LLA_order

c      xbar = 1-x
c      factor=factor*eps*xbar**(eps-1.0)

      if( LLA_order > 0 ) then
        f = 1.0D0 + g1*eps
        
        x_2 = x*x
        if( xxbar > 0) f=f-(1.0D0-x_2)/(2*xxbar)
        if(LLA_order > 1 ) then
          f=f + g2*eps**2
          if( xxbar>0 .and. xbar>0 .and. x>xmin ) then
            log_x=log(x)
            log_xbar = log(xbar)
            f=f - ((1.0D0+3.0D0*x_2)*log_x + 
     >            xbar*( 4.0D0*(1.D0+x)*log_xbar + 5.0D0+x)) 
     >       / (8.0D0*xxbar)*eps
          endif
          if( LLA_order > 2 ) then
            f=f+g3*eps**3
            if( xxbar>0 .and. xbar>0 .and. x>xmin ) then
c            if( xxbar>0 .and. xbar>1110 .and. x>xmin ) then
              f = f - ((1.D0+x) * xbar 
     >          * (6.D0 * Li2_ISR(x) + 12.D0 * log_xbar**2 - 3.D0*pi**2)
     >        + 1.5D0 * (1.0D0 + 8.0D0*x + 3.0D0*x_2) * log_x 
     >        + 6.D0 * (x+5.D0) * xbar * log_xbar 
     >        + 12.D0 * (1.D0+x_2) * log_x * log_xbar 
     >        - (1.D0 + 7.D0*x_2) * log_x**2 / 2.D0 
     >        + (39.D0 - 24.D0*x - 15.D0*x_2) / 4.D0) 
     >        / ( 48.D0 * xxbar) * eps**2
            end if
          end if
        end if
        factor = factor * f
      end if
    
      return 
      end

          
      double precision function Li2_ISR(x)
      double precision x
      double precision pi, pi2_6
      double precision Li2_restricted_ISR
      parameter (pi=3.14159265358979D0)
      parameter (pi2_6=pi**2/6.0D0)

      if (abs(1.0D0-x) < 1.0D-13) then
        Li2_ISR = pi2_6
      else if (abs(1.0D0-x) <  0.5D0) then
        Li2_ISR = pi2_6 - log(1.0D0-x) * log(x) 
     >          - Li2_restricted_ISR (1.0D0-x)
      else if (abs(x).gt.1.d0) then
        print *,'Error in Li2_ISR '
        print *,'Dilogarithm called outside of defined range.'
c       call msg_bug (" Dilogarithm called outside of defined range.")
c       Li2 = -pi2_6 - 0.5_default * log(-x) * log(-x) - Li2_restricted (1/x)
      else
        Li2_ISR = Li2_restricted_ISR (x)
      end if
      
      return 
      end

      double precision function Li2_restricted_ISR (x)
      double precision x
      double precision tmp, z, z2

      z = - log (1-x)
      z2 = z**2
! Horner's rule for the powers z^3 through z^19
      tmp = 43867.0D0/798.0D0
      tmp = tmp * z2 /342.0D0 - 3617.0D0/510.0D0
      tmp = tmp * z2 /272.0D0 + 7.0D0/6.0D0
      tmp = tmp * z2 /210.0D0 - 691.0D0/2730.0D0
      tmp = tmp * z2 /156.0D0 + 5.0D0/66.0D0
      tmp = tmp * z2 /110.0D0 - 1.0D0/30.0D0
      tmp = tmp * z2 / 72.0D0 + 1.0D0/42.0D0
      tmp = tmp * z2 / 42.0D0 - 1.0D0/30.0D0
      tmp = tmp * z2 / 20.0D0 + 1.0D0/6.0D0
! The first three terms of the power series
      Li2_restricted_ISR = z2 * z * tmp / 6.0D0 - 0.25D0 * z2 + z

      return
      end

 
