C* (Update Record)
C*    96/08/02  M. NOJIRI     
C*
      SUBROUTINE GTSFMS3(AMU,TANB,ALF,S2W,AZM,FM3,ASFT,SFM,SFM3,SFM3M,
     &SPHI,IS)
 
      IMPLICIT    REAL*4 ( A-H, O-Z )
      REAL*4      AMU, TANB, ALF, S2W, SFM(7),SFM3(7),ASFT(3),FM3(3)
      REAL*4      C(7), SPHI(3),SFM3M(7)
      INTEGER*4 TYPFRM
      REAL*4    SFMAS(2)
      REAL*4 ATERM
C      ASFT=(ATAU,ATOP,ABOTTOM)
C      SFM(3)=(AMSNL,AMSEL,AMSER,AMSUL,AMSUR,AMSDL,AMSDR)
C      SFM3M=(AMSNL3,AMSTAU1,AMSTAU2,AMST1,AMST2,AMSB1,AMSB2)

      DATA NCALL /0/
      IF(IS.EQ.0) then 
C     calculate SFM3(7) from SFM(7) (extract D term contribution 
C     out) 
C     IF IS.NEQ.0 IT USE THE VALUE AS IT IS INPUTTED.
      C(1) = -0.5
      C(2) =  0.5 -      S2W
      C(3) =             S2W      
      C(4) = -0.5 + 2/3.*S2W
      C(5) =      - 2/3.*S2W
      C(6) =  0.5 - 1/3.*S2W
      C(7) =        1/3.*S2W
      AMZCS = AZM*AZM*(TANB-1)*(TANB+1)/(TANB*TANB+1)
      DO 20 I = 1, 7
         SFM3(I) = SQRT(SFM(I)**2 - C(I)*AMZCS)
         write(*,*) 'sfm3(',I ,sfm3(I)
20    CONTINUE


      ENDIF
      DO I=1,3
      ATERM=ASFT(I)
      SFTMSL=SFM3(2*I)
      SFTMSR=SFM3(2*I+1)
      FRMMAS=FM3(I)   
      write(*,*) aterm,sftmsl,sftmsr,frmmas
      IRET=0 
      WM=azm*sqrt(1-s2w) 
      CALL SFMMIX( WM, S2W,
     &                   AMU, TANB, ATERM,
     &                   SFTMSL, SFTMSR, FRMMAS, I,
     &                   SFMAS, PHIF, THTF,
     &                   IRET  )
      IF(IRET.LT.0) THEN 
        WRITE(6,*) 'UNDEFINED THIRD GENERATION SFMASS'
      ENDIF
      SFM3M(2*I) = SFMAS(1)
      SFM3M(2*I+1) = SFMAS(2)
      SPHI(I) = THTF
      ENDDO
      SFM3M(1)=SQRT(SFM3(1)**2 + C(1)*AMZCS)
      RETURN
      END
