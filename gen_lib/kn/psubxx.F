      SUBROUTINE PSUBxx(P1,P2 , PSUB)                                           
C                                                                               
C This subroutine computes the subtraction of two real vectors.                 
C                                                                               
C INPUT:                                                                        
C       real    P1(0:3)        : real vector                                    
C       real    P2(0:3)        : real vector                                    
C                                                                               
C OUTPUT:                                                                       
C       real    PSUB(0:3)      : subtraction P1-P2                              
C                                                                               
      REAL    P1(0:3),P2(0:3),PSUB(0:3)                                         
C                                                                               
      PSUB(0)=P1(0)-P2(0)                                                       
      PSUB(1)=P1(1)-P2(1)                                                       
      PSUB(2)=P1(2)-P2(2)                                                       
      PSUB(3)=P1(3)-P2(3)                                                       
      RETURN                                                                    
      END                                                                       
