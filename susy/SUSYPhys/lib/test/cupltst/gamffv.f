      subroutine gamffv(amc,amn,amw,gncw,gam)
      
      real*4     amc, amn, amw, gam
      complex*8  gncw(2)
      beta(x1,x2) = sqrt( 1 - 2*(x1+x2) + (x1-x2)**2 )
      
      if ( amc.le.amn+amw ) then
         gam = 0
         return
      endif
      pi  = acos(-1.)
      fct = beta( (amn/amc)**2, (amw/amc)**2 )/16/pi/amc
      tta = ( ( amc**2 + (amn-amw)*(amn+amw) )*amw**2
     .      + ( (amc-amn)*(amc+amn) + amw**2 )*
     .               ( (amc-amn)*(amc+amn) - amw**2 ) )
     .      * ( abs(gncw(1))**2 + abs(gncw(2))**2 )/2/amw**2
     .      - 6*amc*amn*REAL(gncw(1)*CONJG(gncw(2)))
     
      gam = fct*tta
      return
      end
