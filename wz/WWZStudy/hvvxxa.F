      SUBROUTINE HVVXXA(V1,V2,G1,G2,G3,SMASS,SWIDTH , HVV)
C
C This subroutine computes an off-shell scalar current from the vector- 
C vector-scalar coupling.                                               
C                                                                       
C INPUT:                                                                
C       complex V1(6)          : first  vector                        V1
C       complex V2(6)          : second vector                        V2
C       real    G1             : coupling constant                 GVVHa
C       real    G2             : coupling constant                 GVVHb
C       real    G3             : coupling constant                GVVHbt
C       real    SMASS          : mass  of OUTPUT scalar S               
C       real    SWIDTH         : width of OUTPUT scalar S               
C                                                                       
C OUTPUT:                                                               
C       complex HVV(3)         : off-shell scalar current     J(S:V1,V2)
C
      COMPLEX V1(6),V2(6),HVV(3),D,V1V2,P1V2,P2V1
      REAL    Q(0:3),G1,G2,G3,SMASS,SWIDTH,Q2,P1(0:3),P2(0:3)
C
      HVV(2) = V1(5)+V2(5)
      HVV(3) = V1(6)+V2(6)
C
      Q(0)=REAL( HVV(2))
      Q(1)=REAL( HVV(3))
      Q(2)=AIMAG(HVV(3))
      Q(3)=AIMAG(HVV(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      D=-1/CMPLX( Q2-SMASS**2 , MAX(SIGN( SMASS*SWIDTH ,Q2),0.) )
C
      P1P2   = P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)
      P1V2   = P1(0)*V2(1)-P1(1)*V2(2)-P1(2)*V2(3)-P1(3)*V2(4)
      P2V1   = P2(0)*V1(1)-P2(1)*V1(2)-P2(2)*V1(3)-P2(3)*V1(4)
      V1V2   = V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4)
C
      P1(0)=REAL( V1(5))
      P1(1)=REAL( V1(6))
      P1(2)=AIMAG(V1(6))
      P1(3)=AIMAG(V1(5))
C
      P2(0)=REAL( V2(5))
      P2(1)=REAL( V2(6))
      P2(2)=AIMAG(V2(6))
      P2(3)=AIMAG(V2(5))
C
      HVV(1) = G1*V1V2
      HVV(1) = HVV(1)+G2*(P1P2*V1V2-P1V2*P2V1)
      HVV(1) = HVV(1)+G3*(P1(0)*V1(2)*P2(2)*V2(4)
     .                   -P1(0)*V1(2)*P2(3)*V2(3)
     .                   -P1(0)*V1(3)*P2(1)*V2(4)
     .                   +P1(0)*V1(3)*P2(3)*V2(2)
     .                   +P1(0)*V1(4)*P2(1)*V2(3)
     .                   -P1(0)*V1(4)*P2(2)*V2(2)
     .                   -P1(1)*V1(1)*P2(2)*V2(4)
     .                   +P1(1)*V1(1)*P2(3)*V2(3)
     .                   +P1(1)*V1(3)*P2(0)*V2(4)
     .                   -P1(1)*V1(3)*P2(3)*V2(1)
     .                   -P1(1)*V1(4)*P2(0)*V2(3)
     .                   +P1(1)*V1(4)*P2(2)*V2(1)
     .                   +P1(2)*V1(1)*P2(1)*V2(4)
     .                   -P1(2)*V1(1)*P2(3)*V2(2)
     .                   -P1(2)*V1(2)*P2(0)*V2(4)
     .                   +P1(2)*V1(2)*P2(3)*V2(1)
     .                   +P1(2)*V1(4)*P2(0)*V2(2)
     .                   -P1(2)*V1(4)*P2(1)*V2(1)
     .                   -P1(3)*V1(1)*P2(1)*V2(3)
     .                   +P1(3)*V1(1)*P2(2)*V2(2)
     .                   +P1(3)*V1(2)*P2(0)*V2(3)
     .                   -P1(3)*V1(2)*P2(2)*V2(1)
     .                   -P1(3)*V1(3)*P2(0)*V2(2)
     .                   +P1(3)*V1(3)*P2(1)*V2(1))
C
      HVV(1) = D*HVV(1)
C
      RETURN
      END
