! interface subroutine to call user_spectrum 
! prepared for Whizard 

      subroutine whizard_spectrum(isrbm, roots, z1, z2, ebm, &
                           dpdebm, embm, epbm )

      use user, only : spectrum_double

!-- variables for Tim's luminosity function
      real(4) :: roots
      real(8) :: z1,z2
      complex(8), dimension(-2:2,-2:2,2), save :: rhoin_a6f, rhoout_a6f
      real(8) :: dpdebm_dbl 
      real(8), dimension(2) :: sqrts_a6f, x_a6f
      integer, dimension(2), save :: inprt_a6f, mod_a6f
      integer, save :: iflag = 0
      data inprt_a6f/11, -11/
      data rhoin_a6f/50*0.0/
      data rhoout_a6f/50*0.0/

      if ( iflag .eq. 0 ) then
        do 100 i=-2,2
          rhoin_a6f(i,i,1)=1
          rhoin_a6f(i,i,2)=1
100     continue
        iflag = 1
      endif

      mod_a6f(1)=mod(isrbm,100)
      mod_a6f(2)=mod(isrbm,100)
      sqrts_a6f(1)=ROOTS
      sqrts_a6f(2)=ROOTS
      x_a6f(1)=z1
      x_a6f(2)=z2
      
!      print *,'whizard spectrum was called ... isrbm=',isrbm
!      print *,' roots=',roots,' z1,z2=',z1,z2
      call spectrum_double( dpdebm_dbl, x_a6f, inprt_a6f, &
             sqrts_a6f, rhoin_a6f, rhoout_a6f, mod_a6f)
      dpdebm=dpdebm_dbl
      embm=EBM*x_a6f(1)
      epbm=EBM*x_a6f(2)
!      print *,' dpdebm=',dpdebm,' embm, epbm=',embm,epbm

      return
      end subroutine whizard_spectrum

