      use beams, only : strfun_ISR

      implicit real*8(a-h, o-z)
      real*8 factor, x, eps
      integer*4 LLA_order
      logical map

      factor=1.0
      x=0.5

      ene=1000.0
      ame=0.00051
      alphai=137.0
      alpha=1.0/alphai
      pi=acos(-1.0)

      LLA_order=3
      map=.true.

!      do 100 i=0, 10
        i=5
        eps=alpha/pi*2*log(ene/ame)
        x=0.1d0*i
        xx=x
        factor=1
        xfact=1.0
        print *,'=== by isr_function ================================'
        call isr_function(factor, x, eps, LLA_order)
        print *,' i=',i,' x=',x,' factor=',factor

        print *,'==== by strfun_ISR  ==============================='
        call strfun_ISR(xfact, xx, eps, LLA_order, map)
!      print *,' eps=',eps
        print *,' i=',i,' xx=',xx,' xfact=',xfact
100   continue

      print *,' tiny=',tiny(0.2)

      stop
      end

