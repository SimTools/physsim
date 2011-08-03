c ===================================================
      subroutine isr_remnant(x, x0, ran1, ran2, sqrts, pmom)
c
c(Function)
c  generate Pt of ISR photons
c(input Arguments)
c  x : random value, which is returned by strfun_ISR
c  x0 : beam fractional energy before ISR ( x0 after beam strahlung )
c  ran1 : random number to generate pt
c  ran2 : random number to generate photon azimuthal angle
c  sqrts : nominal CM energy in GeV
c(output arguments)
c  pmom(0:3) : Photon momentum.  (E,px,py,pz)
      implicit real*8 (a-h,o-z)
      real*8 pmom(0:3)
      parameter (pi=3.14159265358979D0)
      parameter (twopi=2.0d0*pi)
      parameter (ame=0.000510998902D0)

c ======= Entry Point ==============================

      qmax=sqrts*x0
      xbar=1.0d0-x
      if( xbar .gt. 0.0d0 ) then 
        yran=ran1
        xisr_log=log(1.0d0+(qmax/ame)**2)
        am2=ame**2
        q2 = am2*(1.0d0-xbar*exp(yran*xisr_log))
        tmp = (2.0d0/xbar)*(q2-x*am2)/qmax**2 
      else
        tmp=0.0d0
      endif
      
      if( tmp .ge. 0.0d0 ) then
        st=0.d0
        ct=1.d0
      else
        st2=-tmp*(2.0d0+tmp)
        if( st2 .lt. 0.0d0 ) then
          st=0.d0
          ct=-1.0d0
        else
          st=sqrt(st2)
          ct=1.0d0+tmp
        endif
      endif
      phi=twopi*ran2

cc      call set_beam_remnant(ps, photon, x, i, ct, st, phi)
      E=max((1.0d0-x)*x0*sqrts/2, 0.0d0)
      cphi=cos(phi)
      sphi=sin(phi)
      p=E
      pt=p*st
c      print *,' p=',p,' pt=',pt

      pmom(0)=E
      pmom(1)=pt*cphi
      pmom(2)=pt*sphi
      pmom(3)=p*ct

      return 
      end

