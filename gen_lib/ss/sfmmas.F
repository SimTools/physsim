C=====================================================================C
      SUBROUTINE SFLMAS( COMSCM, G2MAS, TANB, SIN2W, ZM,
     &                   SFMSL, SFMSR, SNEUM,
     &                   SFMSL2, SFMSR2, SNEUM2, IRET )
C=====================================================================C
      IRET=1
C fix sfermion masses
      TAN2B = TANB**2
      COS2B = 1./(1.+TAN2B)
      COSDB = 2.*COS2B-1.
      GMX2  = G2MAS**2 /0.679
      SFMSL2 = COMSCM**2 + 0.52*GMX2 - (0.5-SIN2W)* ZM**2 *COSDB
      SFMSR2 = COMSCM**2 + 0.15*GMX2 - SIN2W* ZM**2 *COSDB
      SNEUM2 = COMSCM**2 + 0.52*GMX2 +  0.5       * ZM**2 *COSDB
 
      IF( SFMSL2.LT.0 .OR. SFMSR2.LT.0 .OR. SNEUM2.LT.0 ) THEN
C       WRITE(6,*) ' G2MAS=',G2MAS, ' SFMSL2=',SFMSL2,
C    &             ' SFMSR2=',SFMSR2, ' SNEUM2=',SNEUM2
        SFMSL2= AMAX1( SFMSL2, 0.)
        SFMSR2= AMAX1( SFMSR2, 0.)
        SNEUM2= AMAX1( SNEUM2, 0.)
        IRET=-1
C       RETURN
      ENDIF
 
      SFMSL = SQRT( SFMSL2 )
      SFMSR = SQRT( SFMSR2 )
      SNEUM = SQRT( SNEUM2 )
 
      RETURN
      END
C=====================================================================C
      SUBROUTINE SFQMAS( COMSCM, G3MAS, TANB, SIN2W, ZM,
     &                   SFDML, SFDMR, SFUML, SFUMR,
     &                   SFDML2, SFDMR2, SFUML2, SFUMR2, IRET )
C=====================================================================C
      IRET=1
C fix sfermion masses
      TAN2B = TANB**2
      COS2B = 1./(1.+TAN2B)
      COSDB = 2.*COS2B-1.
      G3M2  = G3MAS**2
      SFDML2 = COMSCM**2 + 0.837*G3M2 -(0.5-1./3.*SIN2W)*ZM**2 *COSDB
      SFDMR2 = COMSCM**2 + 0.783*G3M2 - 1./3.* SIN2W    *ZM**2 *COSDB
      SFUML2 = COMSCM**2 + 0.837*G3M2 +(0.5-2./3.*SIN2W)*ZM**2 *COSDB
      SFUMR2 = COMSCM**2 + 0.783*G3M2 + 2./3.* SIN2W    *ZM**2 *COSDB
 
      IF( SFDML2.LT.0 .OR. SFDMR2.LT.0   .OR.
     &    SFUML2.LT.0 .OR. SFUMR2.LT.0  ) THEN
C       WRITE(6,*) ' G3MAS=',G3MAS, ' SFDML2=',SFDML2,
C    &             ' SFDMR2=',SFDMR2, ' SFUML2=',SFUML2,
C    &             ' SFUMR2=',SFUMR2
        SFDML2= AMAX1( SFDML2, 0.)
        SFDMR2= AMAX1( SFDMR2, 0.)
        SFUML2= AMAX1( SFUML2, 0.)
        SFUMR2= AMAX1( SFUMR2, 0.)
        IRET=-1
C       RETURN
      ENDIF
 
      SFDML = SQRT( SFDML2 )
      SFDMR = SQRT( SFDMR2 )
      SFUML = SQRT( SFUML2 )
      SFUMR = SQRT( SFUMR2 )
 
      RETURN
      END
