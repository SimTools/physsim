      integer     nxy_max
      parameter ( nxy_max = 50000 )
      real*8  x_L(nxy_max), x_S(nxy_max)
      complex*16  y1_L(nxy_max), y2_L(nxy_max)
      complex*16  y1_S(nxy_max), y2_S(nxy_max)
      integer num_L, num_S
      common/PFIRSTL/ x_L, y1_L, y2_L, num_L
      common/PFIRSTS/ x_S, y1_S, y2_S, num_S
