
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 pmom(0:3)
      integer LLA_order
      logical map

      factor=1.0d0
      sqrts=1000.0d0
      alphai=137.03599967994d0
      alpha=1.0d0/alphai
      pi=3.14159265358979D0
      ame=0.000510998902D0
      LLA_order=3
      map=.true.
      eps=alpha/pi*2*log(sqrts/ame)

      x0=1.0d0
      x=0.1d0
      factor=1.0d0
      call isr_function(factor, x, eps, LLA_order)
      print *,' factor=',factor,' x=',x
   
      call isr_remnant(x, x0, sqrts, pmom)  
      
      print *,' pmom=',(pmom(k),k=0,3)
      
      stop
      end

