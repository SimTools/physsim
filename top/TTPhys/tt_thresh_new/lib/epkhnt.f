CC**********************************************************************
C*
C*======================================================----------===
C*  SUBROUTINE EPKHNT( MODE, ALF,S2W,AMZ,AMW,AMB,
C*                           AMT,ALFS,VTB2,AMH,BTH, EPK,SGPK, NTRY )
C*==================================================--------------===
C*
C* (Purpose)
C*    This subroutine looks for E_peak by scanning through threshold
C*    cross section.
C* (Inputs)
C*       MODE  : (I*4) : mode.
C*       ALF   : (R*8) : alpha(m_Z).
C*       S2W   : (R*8) : sin^2(theta_W).
C*       AMZ   : (R*8) : m_Z.
C*       AMW   : (R*8) : m_W.
C*       AMB   : (R*8) : m_b.
C*       AMT   : (R*8) : m_t.
C*       ALFS  : (R*8) : alpha_s(m_Z).
C*       VTB2  : (R*8) : !V_tb!^2.
C*       AMH   : (R*8) : m_H.
C*       BTH   : (R*8) : beta_H = g_ttH/g_tth(SM).
C*       EPK   : (R*8) : approximate 1S peak position.
C* (Output)
C*       EPK   : (R*8) : Position of 1S peak measured from threshold.
C*       SGPK  : (R*8) : sig_tot at 1S peak.
C*       NTRY  : (I*4) : #trials.
C* (Update Record)
C*    93/04/20  K.Fujii                Original version.
C*
CC**********************************************************************
 
      SUBROUTINE EPKHNT( MODE, ALF,S2W,AMZ,AMW,AMB,
     .                         AMT,ALFS,VTB2,AMH,BTH, EPK,SGPK, NTRY )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       MODE, NTRY
      REAL   *8       AMT, ALFS, VTB2, AMH, BTH, EPK, SGPK
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical and natural constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         GMZ    =  2.5D0
         ESTP1  =  0.010D0
         ESTP2  =  0.001D0
         NTRMAX = 100
      END IF
C--
      VFF   = SQRT(VTB2/1)
C--
C  First decide an approximate peak position.
C--
      EMN = -6.D0
      EMX =  0.D0
      NE  = 30
      DE  = (EMX-EMN)/NE
C--
      MODEP = MODE
      SG0   = 0
      DO 10 IE = 0, NE
         E = EMN + IE*DE
         RS    = 2*AMT + E
         CALL SGTTHR(MODEP,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,
     .                       AMT,VFF,AMH,BTH, RS, SG)
         MODEP = 2
         IF ( SG.GT.SG0 ) THEN
            SG0 = SG
            EPK = E
         ELSE
                                                 GO TO 15
         ENDIF
10    CONTINUE
C--
C  Check which side we are now.
C--
15    RS    = 2*AMT + EPK
      RS    = RS + ESTP2
      CALL SGTTHR(MODEP,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,
     .                       AMT,VFF,AMH,BTH, RS, SG2)
      IF ( SG0.GT.SG2 ) THEN
         DE = -ESTP1
      ELSE
         DE = +ESTP1
      ENDIF
C--
C  Scan through threshold cross section.
C--
      E     = EPK
      SGS   = SG0
      NTRY  = 0
C--
2     NTRY  = NTRY + 1
      E     = E + DE
      RS    = 2*AMT + E
      CALL SGTTHR(MODEP,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,
     .                       AMT,VFF,AMH,BTH, RS, SG)
      IF ( SG.GT.SGS ) THEN
         IF ( NTRY.LT.NTRMAX ) THEN
            SGS   = SG
                                                 GO TO 2
         ELSE
            NTRY  = -1
            RETURN
         ENDIF
      ENDIF
C--
C  E_peak found.
C--
      EL   = E - DE
      EH   = E
      EC   = EL + DE/2
      CALL SGTTHR(MODEP,ALFS,ALF,S2W,AMZ,GMZ,AMW,AMB,
     .                AMT,VFF,AMH,BTH, 2*AMT+EC, SGC )
      SGL  = SGS
      SGH  = SG
      A    = (SGL+SGH-2*SGC)/2
      B    = (SGH-SGL)/2
      XPK  = -B/A/2
      EPK  = EC + XPK*DE/2
      SGPK = SGC - B*B/A/4
C--
C  That's it.
C--
      RETURN
      END
