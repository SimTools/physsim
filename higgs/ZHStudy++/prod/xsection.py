#!/usr/bin/env python
import os
#--
#  Function to calculae xsecton for a given energy range
#--
def subx(emin, emax, n):
    de = (emax-emin)/n
    for i in range(0,n-1):
        e = emin + i*de
        bsfile = "bases." + str(e) + ".root"
        bslog  = "bases." + str(e) + ".log"
        cmd = "jsf -b -q xsection.C'(" + str(e) + ",\"" + bsfile + "\")'" \
		      + " >& " + bslog
        print cmd     
        os.system(cmd)

#--
#  Calculate xsectons
#--
ebd = [198., 230., 270., 530.]
npt = [  16,    8,    13]

for ie in range(0,3):
    emin = ebd[ie]
    emax = ebd[ie+1]
    n    = npt[ie]
    subx(emin, emax, n)

