      SUBROUTINE PSIGNx(P , PSIGN)                                              
C                                                                               
C This subroutine reverses the sign of a four-momentum.                         
C                                                                               
C INPUT:                                                                        
C       real    P(0:3)         : real vector                                    
C                                                                               
C OUTPUT:                                                                       
C       real    PSIGN(0:3)     : reversed four-vector -P                        
C                                                                               
      REAL*8    P(0:3),PSIGN(0:3)                                                 
C                                                                               
      PSIGN(0)=-P(0)                                                            
      PSIGN(1)=-P(1)                                                            
      PSIGN(2)=-P(2)                                                            
      PSIGN(3)=-P(3)                                                            
      RETURN                                                                    
      END                                                                       
