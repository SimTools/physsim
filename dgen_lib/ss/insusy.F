      SUBROUTINE INSUSY(TANB,ALF,S2W,ZM,FM3,
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
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "../ss/include/SSMIX.inc"
      REAL*8    AM0, AMU, AM2, TANB, ALF, ALFS, S2W, ZM, FM3(3),
     .          SPHI(3), ALPH
      COMPLEX*16 GNCW(2,4,2), GCCZ(2,2,2), GNNZ(2,4,4)
      COMPLEX*16 GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     .           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*16 GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     .           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
      COMPLEX*16 GCSNL3(2,2), GCSEL3(2,2), GCSER3(2,2),
     .           GCSDL3(2,2), GCSDR3(2,2), GCSUL3(2,2), GCSUR3(2,2)
      COMPLEX*16 GNSNL3(2,4), GNSEL3(2,4), GNSER3(2,4),
     .           GNSUL3(2,4), GNSUR3(2,4), GNSDL3(2,4), GNSDR3(2,4)
      COMPLEX*16 GCNSE1(2,2), GCNSE2(2,2), 
     .           GCUSD1(2,2), GCUSD2(2,2), GCDSU1(2,2), GCDSU2(2,2)
      COMPLEX*16 GNESE1(2,4), GNESE2(2,4),
     .           GNUSU1(2,4), GNUSU2(2,4), GNDSD1(2,4), GNDSD2(2,4),
     .           GCCH(2,2,2,3), GCNH(2,2,4), GNNH(2,4,4,3) 
C
C==================< Entry Point >======================================
C
C--
C  Set MSSM parameters.
C--
      WM    = ZM*SQRT(1-S2W)
C--
      SM    = AM0
      HM1   = AMU
      G2MAS = AM2
      G3MAS = G2MAS*(ALFS/ALF)*S2W
C--
C  Calculate ino-ino-gauge boson couplings.
C--
      CALL  INOGBS( ALF, S2W,
     .              PHIL, PHIR, ER, ON, ETA,
     .              GNCW, GCCZ, GNNZ )
C--
C  Calculate ino-sfermion-fermion couplings.
C--
C for first and second generations
      CALL  INOFSF( ALF, S2W, WM,
     .              0.D0, 0.D0, 0.D0,
     .              TANB, PHIL, PHIR, ER, ON, ETA,
     .              GCESNL, GCNSEL, GCNSER,
     .              GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     .              GNNSNL, GNESEL, GNESER,
     .              GNUSUL, GNUSUR, GNDSDL, GNDSDR )
C for the third generation
      CALL  INOFSF( ALF, S2W, WM,
     .              FM3(1), FM3(2), FM3(3),
     .              TANB, PHIL, PHIR, ER, ON, ETA,
     .              GCSNL3, GCSEL3, GCSER3,
     .              GCSDL3, GCSDR3, GCSUL3, GCSUR3,
     .              GNSNL3, GNSEL3, GNSER3,
     .              GNSUL3, GNSUR3, GNSDL3, GNSDR3 )
C       write(*,*) 'gnser3'
C       do i=1,4
C       write(*,*) gnser3(1,I),gnser3(2,I)
C       enddo
C with sfermion left-right mixing
      CALL  INOFSF3( SPHI,
     &                    GCSEL3, GCSER3,
     &                    GCSDL3, GCSDR3, GCSUL3, GCSUR3,
     &                    GNSEL3, GNSER3,
     &                    GNSUL3, GNSUR3, GNSDL3, GNSDR3,
     &                    GCNSE1, GCNSE2,
     &                    GCUSD1, GCUSD2, GCDSU1, GCDSU2,
     &                    GNESE1, GNESE2,
     &                    GNUSU1, GNUSU2, GNDSD1, GNDSD2 )
            
C       write(*,*) 'output from inofsf3'
C       do i =1,4
C       write(*,*) gnese1(1,I),gnese1(2,I)
C       enddo
C-- 
C   calculate ino-ino-higgs coupling
C--
      TANA= TAN(ALPH)
      CALL  INOHGS( ALF, S2W, TANB, TANA,
     &                    PHIL, PHIR, ER, ON, ETA,
     &                    GCCH, GCNH, GNNH )
C
C--
C  That's it.
C--
      RETURN
      END












