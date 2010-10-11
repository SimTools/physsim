      subroutine gamsff(ams,am2,am3,g32s,gam)
      
      real*4     ams, am2, am3, gam
      complex*8  g32s(2)
      beta(x1,x2) = sqrt( 1 - 2*(x1+x2) + (x1-x2)**2 )

      if ( ams.le.am2+am3 ) then
         gam = 0
         return
      endif
      pi  = acos(-1.)
      fct = beta( (am2/ams)**2, (am3/ams)**2 )/16/pi/ams
      tta = ( ams**2 - am2**2 -am3**2 )
     .         * ( abs(g32s(1))**2 + abs(g32s(2))**2 )  
     .      -4 * am2*am3 * REAL( g32s(1)*CONJG(g32s(2)) )
      gam = fct*tta
      return
      end
