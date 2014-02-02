#!/usr/bin/env python
import os
#--
#  Function to calculae xsecton for a given energy range
#--
def subx(emin, emax, n):
    de = (emax-emin)/n
    for i in range(0,n):
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
#ebd = [320., 350., 400., 1200.]
ebd = [220., 250., 300., 1100.]
npt = [  15,   10,  16]

for ie in range(0,3):
    emin = ebd[ie]
    emax = ebd[ie+1]
    n    = npt[ie]
    subx(emin, emax, n)

