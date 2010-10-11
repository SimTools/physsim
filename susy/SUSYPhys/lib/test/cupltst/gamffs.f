      subroutine gamffs(am1,am2,ams,g21s,gam)
      
      real*4     am1, am2, ams, gam
      complex*8  g21s(2)
      beta(x1,x2) = sqrt( 1 - 2*(x1+x2) + (x1-x2)**2 )

      if ( am1.le.am2+ams ) then
         gam = 0
         return
      endif
      pi  = acos(-1.)
      fct = beta( (am2/am1)**2, (ams/am1)**2 )/16/pi/am1
      tta = ( am1**2 + am2**2 -ams**2 )
     .         * ( abs(g21s(1))**2 + abs(g21s(2))**2 )/2  
     .      + 2 * am1*am2 * REAL( g21s(1)*CONJG(g21s(2)) )
      gam = fct*tta
      return
      end
