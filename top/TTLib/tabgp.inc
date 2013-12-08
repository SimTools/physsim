#include "nmax.inc"
      integer     nq_max
      parameter ( nq_max = 50000 )
      complex*8  array_GP(0:n_max,0:nq_max)
      common/tabulateGP/array_GP
