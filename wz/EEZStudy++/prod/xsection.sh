#!/bin/sh
for ((npt=1; $npt != 4; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((90+5*$npt))\)
done
for ((npt=0; $npt != 8; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((110+20*$npt))\)
done
for ((npt=0; $npt != 13; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((270+100*$npt))\)
done
