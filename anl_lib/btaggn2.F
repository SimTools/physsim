      SUBROUTINE BTAGGN(NTRK,IQJT,IPTOJ,IPTDEX,IGNDEX,
     .                  LNxDT,RDATA,GDATA,NBTKCT,ISPT)
 
      INTEGER*4  NTRK, IQJT, IPTOJ(*), IPTDEX(*), IGNDEX(*),
     .           LNxDT, NBTKCT, ISPT
      REAL   *4  RDATA(LNxDT,*), GDATA(LNxDT,*)
C--
C  Tagging efficiency for 2.5-sigma cut.
C     EFF(NTK,1) = tagging efficiency for charm.
C        (NTK,2) = tagging efficiency for bottom.
C--
      REAL   *4  EFF(0:3,2)
      DATA EFF   / 1.00, 0.70, 0.38, 0.14,
     .             1.00, 0.90, 0.78, 0.63 /
C--
C  Fudge factor to make the resultant tagging efficiency consistent
C  with the input.
C--
      DATA FCR   / 1.24 /
C--
      DATA ISEED / 3141 /
C
C==================<<  Entry Point  >>================================
C
C--
C  First decide origin of this jet.
C--
      CALL QTAGGN(NTRK,IQJT,IPTOJ,IPTDEX,IGNDEX,
     .            LNxDT,RDATA,GDATA,NBTKCT,ISPT)
      IF ( ISPT.EQ.0 )                           GO TO 9999
C--
C  This jet is heavy-flavored. Reset ISPT according to tagging
C  efficiencies.
C--
      ICT = MIN(MAX(NBTKCT,0),3)
      ICB = IABS(ISPT)
      IF ( RAN(ISEED).GT.EFF(ICT,ICB)*FCR ) ISPT = 0
C--
C  That's it.
C--
9999  RETURN
      END
 
      SUBROUTINE QTAGGN(NTRK,IQJT,IPTOJ,IPTDEX,IGNDEX,
     .                  LNxDT,RDATA,GDATA,NBTKCT,ISPT)
 
      INTEGER*4  NTRK, IQJT, IPTOJ(*), IPTDEX(*), IGNDEX(*),
     .           LNxDT, NBTKCT, ISPT
      REAL   *4  RDATA(LNxDT,*), GDATA(LNxDT,*)
C--
C  Weak-decaying c-hadrons.
C--
      PARAMETER ( NCHD = 6 )
      INTEGER*4  IDCHD(NCHD)
      DATA IDCHD /    411,    421,    431,
     .               4122,   4322,   4312 /
C--
C  Weak-decaying b-hadrons.
C--
      PARAMETER ( NBHD = 17 )
      INTEGER*4  IDBHD(NBHD)
      DATA IDBHD /    511,    521,    531,
     .             100145, 100146, 100147, 100148, 100149, 100150,
     .             100241, 100242, 100243,
     .             100293, 100294, 100295, 100296, 100297 /
C
C==================<<  Entry Point  >>================================
C
C--
C  Reset return flag.
C--
      ISPT = 0
C--
C  Loop over charged tracks in jet #IQJT.
C--
      NBTK  = 0
      NCTK  = 0
      EBMX  = 0
      ECMX  = 0
      DO 2000 ITK = 1, NTRK
         IJT = IPTOJ(ITK)
         IF ( IJT.NE.IQJT )                      GO TO 2000
         IPT = IPTDEX(ITK)
         ICH = RDATA(9,IPT)
         IF ( ICH.EQ.0 )                         GO TO 2000
         ETK   = RDATA(4,IPT)
         IGN   = RDATA(8,IPT)
C--
C  Find parent b-hadrons.
C--
10       IPR = GDATA(14,IGNDEX(IGN))
C--
C  Reached a parent parton.
C--
         IF ( IPR.LT.0 ) THEN
                                                 GO TO 2000
C--
C  It's not a parton yet. Check if it's a c/b-hadron.
C--
         ELSE
            ID  = GDATA(2,IGNDEX(IPR))
            IDA = IABS(ID)
C--
            DO 100 ICHD = 1, NCHD
               IF ( IDCHD(ICHD).EQ.IDA ) THEN
                  NCTK = NCTK + 1
                  IF ( ETK.GT.ECMX ) THEN
                     ECMX = ETK
                     IF ( ID.GT.0 ) THEN
                        ISPT = +1
                     ELSE
                        ISPT = -1
                     ENDIF
                  ENDIF
                                                 GO TO 250
               ENDIF
100         CONTINUE
C--
            DO 200 IBHD = 1, NBHD
               IF ( IDBHD(IBHD).EQ.IDA ) THEN
                  NBTK = NBTK + 1
                  IF ( ETK.GT.EBMX ) THEN
                     EBMX = ETK
                     IF ( ID.GT.0 ) THEN
                        ISPT = +2
                     ELSE
                        ISPT = -2
                     ENDIF
                  ENDIF
                                                 GO TO 2000
               ENDIF
200         CONTINUE
C--
C  This is not a b-hadron. Try its parent.
C--
250         IGN = IPR
                                                 GO TO 10
         ENDIF
2000  CONTINUE
C--
C  Check if enough # of b-originated tracks there.
C--
      IF ( NBTK.LT.NBTKCT .AND. NCTK.LT.NBTKCT ) ISPT = 0
C--
C  That's it.
C--
9999  RETURN
      END
