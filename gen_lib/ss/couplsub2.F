      SUBROUTINE  INOFSF3( SPHI,
     &                    GCNSEL, GCNSER,
     &                    GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     &                    GNESEL, GNESER,
     &                    GNUSUL, GNUSUR, GNDSDL, GNDSDR,
     &                    GCNSE1, GCNSE2,
     &                    GCUSD1, GCUSD2, GCDSU1, GCDSU2,
     &                    GNESE1, GNESE2,
     &                    GNUSU1, GNUSU2, GNDSD1, GNDSD2 )
C    GCNSEL(i,j) (i=1..2, j=1..2 ) / chargino(j)-neutrino-Selectron_L
C    GCNSER(i,j) (i=1..2, j=1..2 ) / chargino(j)-neutrino-Selectron_R
C    GCDSUL(i,j) (i=1..2, j=1..2 ) / chargino(j)-d_quark-Su_quark_L
C    GCDSUR(i,j) (i=1..2, j=1..2 ) / chargino(j)-d_quark-Su_quark_R
C    GCUSDL(i,j) (i=1..2, j=1..2 ) / chargino(j)-u_quark-Sd_quark_L
C    GCUSDR(i,j) (i=1..2, j=1..2 ) / chargino(j)-u_quark-Sd_quark_R
C    GNNSNL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-neutrino-Sneutrino_L
C    GNESEL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-electron-Selectron_L
C    GNESER(i,j) (i=1..2, j=1..4 ) / neutralino(j)-electron-Selectron_R
C    GNUSUL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_L
C    GNUSUR(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_R
C    GNDSDL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_L
C    GNDSDR(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_R
C              i=1/left-handed-coupling,  i=2/right-handed
C-------------------------------------------------------------------
C    GCNSE1(i,j) (i=1..2, j=1..2 ) / chargino(j)-neutrino-Selectron_1
C    GCNSE2(i,j) (i=1..2, j=1..2 ) / chargino(j)-neutrino-Selectron_2
C    GCDSU1(i,j) (i=1..2, j=1..2 ) / chargino(j)-d_quark-Su_quark_1
C    GCDSU2(i,j) (i=1..2, j=1..2 ) / chargino(j)-d_quark-Su_quark_2
C    GCUSD1(i,j) (i=1..2, j=1..2 ) / chargino(j)-u_quark-Sd_quark_1
C    GCUSD2(i,j) (i=1..2, j=1..2 ) / chargino(j)-u_quark-Sd_quark_2
C    GNESE1(i,j) (i=1..2, j=1..4 ) / neutralino(j)-electron-Selectron_1
C    GNESE2(i,j) (i=1..2, j=1..4 ) / neutralino(j)-electron-Selectron_2
C    GNUSU1(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_1
C    GNUSU2(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_2
C    GNDSD1(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_1
C    GNDSD2(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_2
C-----------------------
C  BASIC RELATION
C  G1 =  GL COS(THETA) + GR SIN(THETA)
C  G2 = -GL SIN(THETA) + GR COS(THETA)
C----------------------------------------------------------------------C
      COMPLEX*8  GCNSEL(2,2), GCNSER(2,2),
     &          GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*8  GNESEL(2,4), GNESER(2,4),
     &          GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
      COMPLEX*8 GCNSE1(2,2), GCNSE2(2,2),
     &          GCUSD1(2,2), GCUSD2(2,2), GCDSU1(2,2), GCDSU2(2,2)
      COMPLEX*8 GNESE1(2,4), GNESE2(2,4),
     &          GNUSU1(2,4), GNUSU2(2,4), GNDSD1(2,4), GNDSD2(2,4)
      REAL*4 SPHI(3)
C--
C     INITIALIZE COS OR SIN (THETA_SF)
C--      
      SSE = SIN(SPHI(1))
      CSE = COS(SPHI(1))
      SSU = SIN(SPHI(2))
      CSU = COS(SPHI(2))
      SSD = SIN(SPHI(3))
      CSD = COS(SPHI(3))
C--
      do i = 1,2
      do j = 1,2
C  BASIC RELATION
C  G1 =  GL COS(THETA) + GR SIN(THETA)
C  G2 = -GL SIN(THETA) + GR COS(THETA)
        GCNSE1(i,j) =  CSE * GCNSEL(I,J) + SSE * GCNSER(I,J)
        GCNSE2(i,j) = -SSE * GCNSEL(I,J) + CSE * GCNSER(I,J)
        GCDSU1(i,j) =  CSU * GCDSUL(I,J) + SSU * GCDSUR(I,J)
        GCDSU2(i,j) = -SSU * GCDSUL(I,J) + CSU * GCDSUR(I,J)
        GCUSD1(i,j) =  CSD * GCUSDL(I,J) + SSD * GCUSDR(I,J)
        GCUSD2(i,j) = -SSD * GCUSDL(I,J) + CSD * GCUSDR(I,J)
      ENDDO
      ENDDO
      DO I=1,2
      DO J=1,4
        GNESE1(i,j) =  CSE * GNESEL(I,J) + SSE * GNESER(I,J)
        GNESE2(i,j) = -SSE * GNESEL(I,J) + CSE * GNESER(I,J)
        GNUSU1(i,j) =  CSU * GNUSUL(I,J) + SSU * GNUSUR(I,J)
        GNUSU2(i,j) = -SSU * GNUSUL(I,J) + CSU * GNUSUR(I,J)
        GNDSD1(i,j) =  CSD * GNDSDL(I,J) + SSD * GNDSDR(I,J)
        GNDSD2(i,j) = -SSD * GNDSDL(I,J) + CSD * GNDSDR(I,J)
      ENDDO
      ENDDO
      END
C-----------------------
        
    







