#!/bin/sh
for ((npt=0; $npt != 6; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((270+2*$npt))\)
done
for ((npt=0; $npt != 23; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((280+10*$npt))\)
done
