C   22/05/92 205220052  MEMBER NAME  FVIXXX   (S)           FORTRAN
C
C ----------------------------------------------------------------------
C
      SUBROUTINE FVIXXC(FI,VC,G,FMASS,FWIDTH , FVI)
C
C This subroutine computes an off-shell fermion wavefunction from a
C flowing-IN external fermion and a vector boson.
C
C INPUT:
C       complex FI(6)          : flow-in  fermion                   !FI>
C       complex VC(6)          : input    vector                      V
C       complex G(2)           : coupling constants                  GVF
C       real    FMASS          : mass  of OUTPUT fermion F'
C       real    FWIDTH         : width of OUTPUT fermion F'
C
C OUTPUT:
C       complex FVI(6)         : off-shell fermion             !F',V,FI>
C
      COMPLEX FI(6),VC(6),FVI(6),SL1,SL2,SR1,SR2,D,CI,G(2)
      REAL    PF(0:3),FMASS,FWIDTH,PF2
C
      FVI(5) = FI(5)-VC(5)
      FVI(6) = FI(6)-VC(6)
C
      PF(0)=REAL( FVI(5))
      PF(1)=REAL( FVI(6))
      PF(2)=AIMAG(FVI(6))
      PF(3)=AIMAG(FVI(5))
      CI=CMPLX(0.,1.)
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      D=-1./CMPLX( PF2-FMASS**2 , MAX(SIGN( FMASS*FWIDTH ,PF2),0.) )
      SL1= (VC(1)+   VC(4))*FI(1)
     &    +(VC(2)-CI*VC(3))*FI(2)
      SL2= (VC(2)+CI*VC(3))*FI(1)
     &    +(VC(1)-   VC(4))*FI(2)
C
      IF (ABS(G(2)).EQ.0.) GOTO 10
C
      SR1= (VC(1)-   VC(4))*FI(3)
     &    -(VC(2)-CI*VC(3))*FI(4)
      SR2=-(VC(2)+CI*VC(3))*FI(3)
     &    +(VC(1)+   VC(4))*FI(4)
C
      FVI(1) = ( G(1)*((PF(0)-PF(3))*SL1 -CONJG(FVI(6))*SL2)
     &          +G(2)*FMASS*SR1)*D
      FVI(2) = ( G(1)*(      -FVI(6)*SL1 +(PF(0)+PF(3))*SL2)
     &          +G(2)*FMASS*SR2)*D
      FVI(3) = ( G(2)*((PF(0)+PF(3))*SR1 +CONJG(FVI(6))*SR2)
     &          +G(1)*FMASS*SL1)*D
      FVI(4) = ( G(2)*(       FVI(6)*SR1 +(PF(0)-PF(3))*SR2)
     &          +G(1)*FMASS*SL2)*D
C
      RETURN
C
  10  FVI(1) = G(1)*((PF(0)-PF(3))*SL1 -CONJG(FVI(6))*SL2)*D
      FVI(2) = G(1)*(      -FVI(6)*SL1 +(PF(0)+PF(3))*SL2)*D
      FVI(3) = G(1)*FMASS*SL1*D
      FVI(4) = G(1)*FMASS*SL2*D
C
      RETURN
      END
C   22/05/92 205220052  MEMBER NAME  FVOXXX   (S)           FORTRAN
C
C ----------------------------------------------------------------------
C
      SUBROUTINE FVOXXC(FO,VC,G,FMASS,FWIDTH , FVO)
C
C This subroutine computes an off-shell fermion wavefunction from a
C flowing-OUT external fermion and a vector boson.
C
C INPUT:
C       complex FO(6)          : flow-out fermion                   <FO!
C       complex VC(6)          : input    vector                      V
C       complex G(2)           : coupling constants                  GVF
C       real    FMASS          : mass  of OUTPUT fermion F'
C       real    FWIDTH         : width of OUTPUT fermion F'
C
C OUTPUT:
C       complex FVO(6)         : off-shell fermion             <FO,V,F'!
C
      COMPLEX FO(6),VC(6),FVO(6),SL1,SL2,SR1,SR2,D,CI,G(2)
      REAL    PF(0:3),FMASS,FWIDTH,PF2
C
      FVO(5) = FO(5)+VC(5)
      FVO(6) = FO(6)+VC(6)
C
      PF(0)=REAL( FVO(5))
      PF(1)=REAL( FVO(6))
      PF(2)=AIMAG(FVO(6))
      PF(3)=AIMAG(FVO(5))
      CI=CMPLX(0.,1.)
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      D=-1./CMPLX( PF2-FMASS**2 , MAX(SIGN( FMASS*FWIDTH ,PF2),0.) )
      SL1= (VC(1)+   VC(4))*FO(3)
     &    +(VC(2)+CI*VC(3))*FO(4)
      SL2= (VC(2)-CI*VC(3))*FO(3)
     &    +(VC(1)-   VC(4))*FO(4)
C
      IF (ABS(G(2)).EQ.0.) GOTO 10
C
      SR1= (VC(1)-   VC(4))*FO(1)
     &    -(VC(2)+CI*VC(3))*FO(2)
      SR2=-(VC(2)-CI*VC(3))*FO(1)
     &    +(VC(1)+   VC(4))*FO(2)
C
      FVO(1) = ( G(2)*( (PF(0)+PF(3))*SR1        +FVO(6)*SR2)
     &          +G(1)*FMASS*SL1)*D
      FVO(2) = ( G(2)*( CONJG(FVO(6))*SR1 +(PF(0)-PF(3))*SR2)
     &          +G(1)*FMASS*SL2)*D
      FVO(3) = ( G(1)*( (PF(0)-PF(3))*SL1        -FVO(6)*SL2)
     &          +G(2)*FMASS*SR1)*D
      FVO(4) = ( G(1)*(-CONJG(FVO(6))*SL1 +(PF(0)+PF(3))*SL2)
     &          +G(2)*FMASS*SR2)*D
C
      RETURN
C
  10  FVO(1) = G(1)*FMASS*SL1*D
      FVO(2) = G(1)*FMASS*SL2*D
      FVO(3) = G(1)*( (PF(0)-PF(3))*SL1        -FVO(6)*SL2)*D
      FVO(4) = G(1)*(-CONJG(FVO(6))*SL1 +(PF(0)+PF(3))*SL2)*D
C
      RETURN
      END

      SUBROUTINE VSSXXC(VC,S1,S2,GC , VERTEX)
C
C This subroutine computes an amplitude from the vector-scalar-scalar   
C coupling.  The coupling is absent in the minimal SM in unitary gauge. 
C                                                                       
C       complex VC(6)          : input  vector                        V 
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       complex GC             : coupling constant (S1 charge)          
C                                                                       
C Examples of the coupling constant GC for SUSY particles are as follows:
C   -----------------------------------------------------------         
C   |    S1    | (Q,I3) of S1  ||   V=A   |   V=Z   |   V=W   |         
C   -----------------------------------------------------------         
C   | nu~_L    | (  0  , +1/2) ||   ---   |  GZN(1) |  GWF(1) |         
C   | e~_L     | ( -1  , -1/2) ||  GAL(1) |  GZL(1) |  GWF(1) |         
C   | u~_L     | (+2/3 , +1/2) ||  GAU(1) |  GZU(1) |  GWF(1) |         
C   | d~_L     | (-1/3 , -1/2) ||  GAD(1) |  GZD(1) |  GWF(1) |         
C   -----------------------------------------------------------         
C   | e~_R-bar | ( +1  ,  0  ) || -GAL(2) | -GZL(2) | -GWF(2) |         
C   | u~_R-bar | (-2/3 ,  0  ) || -GAU(2) | -GZU(2) | -GWF(2) |         
C   | d~_R-bar | (+1/3 ,  0  ) || -GAD(2) | -GZD(2) | -GWF(2) |         
C   -----------------------------------------------------------         
C where the S1 charge is defined by the flowing-OUT quantum number.     
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude                Gamma(V,S1,S2)
C
      COMPLEX GC,VC(6),S1(3),S2(3),VERTEX
      REAL    P(0:3)
C
      P(0)=REAL( S1(2)-S2(2))
      P(1)=REAL( S1(3)-S2(3))
      P(2)=AIMAG(S1(3)-S2(3))
      P(3)=AIMAG(S1(2)-S2(2))
C
      VERTEX = GC*S1(1)*S2(1)
     &        *(VC(1)*P(0)-VC(2)*P(1)-VC(3)*P(2)-VC(4)*P(3))
C
      RETURN
      END

      SUBROUTINE IOVXXC(FI,FO,VC,G , VERTEX)
C
C This subroutine computes an amplitude of the fermion-fermion-vector   
C coupling.                                                             
C                                                                       
C INPUT:                                                                
C       complex FI(6)          : flow-in  fermion                   |FI>
C       complex FO(6)          : flow-out fermion                   <FO|
C       complex VC(6)          : input    vector                      V 
C       real    G(2)           : coupling constants                  GVF
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude                     <FO|V|FI>
C
      COMPLEX FI(6),FO(6),VC(6),VERTEX
C>>>
C     REAL    G(2)
      COMPLEX G(2)
C>>>
C
      VERTEX =  G(1)*( (FO(3)*FI(1)+FO(4)*FI(2))*VC(1)
     &                +(FO(3)*FI(2)+FO(4)*FI(1))*VC(2)
     &                -(FO(3)*FI(2)-FO(4)*FI(1))*VC(3)*CMPLX(0.,1.)
     &                +(FO(3)*FI(1)-FO(4)*FI(2))*VC(4)             )
      IF (ABS(G(2)).NE.0.) VERTEX = VERTEX
     &        + G(2)*( (FO(1)*FI(3)+FO(2)*FI(4))*VC(1)
     &                -(FO(1)*FI(4)+FO(2)*FI(3))*VC(2)
     &                +(FO(1)*FI(4)-FO(2)*FI(3))*VC(3)*CMPLX(0.,1.)
     &                -(FO(1)*FI(3)-FO(2)*FI(4))*VC(4)             )
C
      RETURN
      END
