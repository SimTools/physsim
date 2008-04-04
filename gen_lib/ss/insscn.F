      SUBROUTINE INSSCN(ALF,ALFS,WM,ZM,AM0,AMU,AM2,TANB,AMA,ASFT)

      IMPLICIT    REAL*4 ( A-H, O-Z )
#include "../sm/include/SMPTAB.inc"     
#include "../ss/include/SSCONS.inc"
#include "../ss/include/SSPTAB.inc"
#include "../ss/include/SSCUPL.inc"
#include "../ss/include/SSMIX.inc"
      REAL*4   ALF, ALFS, WM, ZM, AM0, AMU, AM2, TANB, AMA
      REAL*4   FM(3), FM3P(3),FM3(3), ASFT(3), xSFM3(7), xSFM(7)
C
C========< Entry Point >================================================
C version without using user defined sfermion masses.
      IS=0
C--
C  Initialize SUSY parameters.
C--
      xAM0     = AM0
      xAMU     = AMU
      xAM2     = AM2
      xTANB    = TANB
      xAMA     = AMA
      xASFT(1) = ASFT(1)
      xASFT(2) = ASFT(2)
      xASFT(3) = ASFT(3)
C--
C  Initialize particle table and coupling constants.
C--
      S2W   = ( 1 - WM/ZM )*( 1 + WM/ZM )
      FM(1) = 0
      FM(2) = 0
      FM(3) = 0

      FM3P(1) = AMF(3,2,1)
      FM3P(2) = AMF(3,1,2)
      FM3P(3) = AMF(3,2,2)
C      CALL RUNMASS(FM3P,FM3,ZM)
      FM3(1)  = FM3P(1)
      FM3(2)  = FM3P(2)
      FM3(3)  = FM3P(3)
C--
      PRINT *, '******** INSSCN ********************************** '
      PRINT *, '    MTAU   = ', FM3P(1)
      PRINT *, '    MTOP   = ', FM3P(2)
      PRINT *, '    MB     = ', FM3P(3)
C      PRINT *, '    mtau   = ', FM3(1)
C      PRINT *, '    mtop   = ', FM3(2)
C      PRINT *, '    mb(mz) = ', FM3(3)
      PRINT *, ' '
      PRINT *, '    AM0    = ', AM0
      PRINT *, '    AMU    = ', AMU
      PRINT *, '    AM2    = ', AM2
      PRINT *, '    TNB    = ', TANB
      PRINT *, '    A_tau  = ', ASFT(1)
      PRINT *, '    A_t    = ', ASFT(2)
      PRINT *, '    A_b    = ', ASFT(3)
      PRINT *, ' '
C--       
      CALL INMASS(ALF,ALFS,WM,ZM,AM0,AMU,AM1,AM2,AM3,TANB,AMA,IS,
     &           xSFM, xSFM3,xASFT,FM3,
     &           SWM,SZM,SFM,SHM,ALPH,SFM3M,SPHI)
C--
C  Print out modified parameters.
C--
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' SUSY Higgs masses                 '
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    Mh     = ', SHM(1)
      PRINT *, '    MH     = ', SHM(2)
      print *, '    MA     = ', SHM(3)
      print *, '    MH+    = ', SHM(4)
C      sinalp= sin(alph)
C      cosalp= cos(alph)
C      print *, 'sinalp,cosalp', sinalp, cosalp 
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' 1st/2nd generation sfermion masses'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    AMSNL  = ', SFM(1)
      PRINT *, '    AMSEL  = ', SFM(2)
      PRINT *, '    AMSER  = ', SFM(3)
      PRINT *, '    AMSUL  = ', SFM(4)
      PRINT *, '    AMSUR  = ', SFM(5)
      PRINT *, '    AMSDL  = ', SFM(6)
      PRINT *, '    AMSDR  = ', SFM(7)
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' 3rd generation sfermion masses'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    AMSNT  = ', SFM3M(1)
      PRINT *, '    AMSE1  = ', SFM3M(2)
      PRINT *, '    AMSE2  = ', SFM3M(3)
      PRINT *, '    AMSU1  = ', SFM3M(4)
      PRINT *, '    AMSU2  = ', SFM3M(5)
      PRINT *, '    AMSD1  = ', SFM3M(6)
      PRINT *, '    AMSD2  = ', SFM3M(7)
      PRINT *, '        '
      PRINT *, ' ----------------------------------'
      PRINT *, ' Chargino masses'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    AMSW1  = ', SWM(1)
      PRINT *, '    AMSW2  = ', SWM(2)
      PRINT *, '        '
      PRINT *, ' ----------------------------------'
      PRINT *, ' Neutralino masses'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    AMSZ1  = ', SZM(1)
      PRINT *, '    AMSZ2  = ', SZM(2)
      PRINT *, '    AMSZ3  = ', SZM(3)
      PRINT *, '    AMSZ4  = ', SZM(4)
      PRINT *, '        '
c
c Add on 05/28/97
c
        print*,' ' 
        print*,'-------  SUSY parameters ------------------------'
        print*,'  M0   =',AM0
        print*,'  M1   =',AM1
        print*,'  M2   =',AM2
        print*,'  mu   =',AMU
        print*,'  tanB =',TANB
	print*,'-------  Mixing angle of third generation -------'
	print*,'   stau mixing    =',sin(sphi(1))
	print*,'   stop mixing    =',sin(sphi(2))
	print*,'   sbottom mixing =',sin(sphi(3))
	print*,'-------------------------------------------------'
        print*,' ' 
C     
         CALL INSUSY(TANB,ALF,S2W,ZM,FM3,
     .                  SPHI, ALPH,
     .                  GNCW,GCCZ,GNNZ,
     .                  GCESNL,GCNSEL,GCNSER,
     .                  GCUSDL,GCUSDR,GCDSUL,GCDSUR,
     .                  GNNSNL,GNESEL,GNESER,
     .                  GNUSUL,GNUSUR,GNDSDL,GNDSDR, 
     .                  GCSNL3,GCNSE1,GCNSE2,
     .                  GCUSD1,GCUSD2,GCDSU1,GCDSU2,
     .                  GNSNL3,GNESE1,GNESE2,
     .                  GNUSU1,GNUSU2,GNDSD1,GNDSD2,
     .                  GCCH,GCNH,GNNH)
C--
       print *, ' '
       print *, '----- Coupling Constants -------'
       print *, ' '
       WRITE(*,*) 'C1-E-NL ', (GCESNL(I,1),I=1,2) 
       WRITE(*,*) 'C1-N-EL ', (GCNSEL(I,1),I=1,2)
       WRITE(*,*) 'C1-N-ER ', (GCNSER(I,1),I=1,2)
       WRITE(*,*) 'C1-E-NL3',(GCSNL3(I,1),I=1,2)
       WRITE(*,*) 'C1-N-E1 ', (GCNSE1(I,1),I=1,2)
       WRITE(*,*) 'C1-N-E2 ', (GCNSE2(I,1),I=1,2)
       WRITE(*,*) 'N0-N-NL ',  (GNNSNL(I,1), I=1,2)
       WRITE(*,*) 'N0-E-EL ',  (GNESEL(I,1), I=1,2)
       WRITE(*,*) 'N0-E-ER ',  (GNESER(I,1), I=1,2)
       WRITE(*,*) 'N0-N-NL3',  (GNSNL3(I,1), I=1,2)
       WRITE(*,*) 'N0-N-E1 ',  (GNESE1(I,1), I=1,2)
       WRITE(*,*) 'N0-N-E2 ',  (GNESE2(I,1), I=1,2)
       print *, ' '
C--
       CALL INHGCN(ZM,WM,TANB,ALPH, SHM, GMSH)
C--
       PRINT *, ' '
      PRINT *, ' ----------------------------------'
       PRINT *, ' SUSY Higgs widths'
      PRINT *, ' ----------------------------------'
       PRINT *, ' '
       PRINT *, '    GMh    = ', GMSH(1)
       PRINT *, '    GMH    = ', GMSH(2)
       PRINT *, '    GMA    = ', GMSH(3)
       PRINT *, '    GMH+   = ', GMSH(4)
       PRINT *, ' '
c
C    print out ino mixing angle for check
C    looks OK. 
C      PRINT '(f6.4,f6.4)', OL(1,1),OL(1,2)
C      PRINT '(F6.4,F6.4)', OL(2,1),OL(2,2)
C      PRINT '(F6.4,F6.4)', OR(1,1),OR(1,2)
C      PRINT '(F6.4,F6.4)', OR(2,1),OR(2,2)
C      PRINT '(F6.4,F6.4,F6.4,F6.4)', ON(1,1),ON(1,2),ON(1,3),ON(1,4)
C      PRINT '(F6.4,F6.4,F6.4,F6.4)', ON(2,1),ON(2,2),ON(2,3),ON(2,4)
C      PRINT '(F6.4,F6.4,F6.4,F6.4)', ON(3,1),ON(3,2),ON(3,3),ON(3,4)
C      PRINT '(F6.4,F6.4,F6.4,F6.4)', ON(4,1),ON(4,2),ON(4,3),ON(4,4)
C      PRINT '(F6.4,F6.4,F6.4,F6.4)', ETA(1),ETA(2),ETA(3),ETA(4)
C--
C      write(*,*) 'input to sf2br'
C      do i=1,4
C      write(*,*) gnese1(1,i),gnese1(2,i)
C      enddo
      CALL SF2BR(  ALF, S2W, WM, ZM, SWM, SZM, SFM, SFM3M, SPHI, 
     &   SHM, FM3P,
     &   GNCW, GCCZ, GNNZ,GCNH,GNNH,
     &   GCESNL, GCNSEL, GCNSER,
     &   GNNSNL, GNESEL, GNESER,
     &   GCSNL3,GCNSE1,GCNSE2,
     &   GNSNL3,GNESE1,GNESE2,
     &   GMSF,GMSW,GMSZ,GMSF3)
C--
C  Print out modified parameters.
C--       
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' 1st/2nd generation sfermion widths'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    GMSNL  = ', GMSF(1)
      PRINT *, '    GMSEL  = ', GMSF(2)
      PRINT *, '    GMSER  = ', GMSF(3)
      PRINT *, '    GMSUL  = ', GMSF(4)
      PRINT *, '    GMSUR  = ', GMSF(5)
      PRINT *, '    GMSDL  = ', GMSF(6)
      PRINT *, '    GMSDR  = ', GMSF(7)
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' 3rd generation sfermion widths'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    GMSNL  = ', GMSF3(1)
      PRINT *, '    GMSE1  = ', GMSF3(2)
      PRINT *, '    GMSE2  = ', GMSF3(3)
      PRINT *, '    GMSU1  = ', GMSF3(4)
      PRINT *, '    GMSU2  = ', GMSF3(5)
      PRINT *, '    GMSD1  = ', GMSF3(6)
      PRINT *, '    GMSD2  = ', GMSF3(7)
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' Chargino widths'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    GMSW1  = ', GMSW(1)
      PRINT *, '    GMSW2  = ', GMSW(2)
      PRINT *, ' '
      PRINT *, ' ----------------------------------'
      PRINT *, ' Neutralino widths'
      PRINT *, ' ----------------------------------'
      PRINT *, ' '
      PRINT *, '    GMSZ1  = ', GMSZ(1)
      PRINT *, '    GMSZ2  = ', GMSZ(2)
      PRINT *, '    GMSZ3  = ', GMSZ(3)
      PRINT *, '    GMSZ4  = ', GMSZ(4)
      PRINT *, '        '
      PRINT *, ' *****************************************************'
C--
C  That's it.
C--
      RETURN
      END
