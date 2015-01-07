#!/bin/sh
for ((npt=0; $npt != 4; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((310+10*$npt))\)
done
for ((npt=0; $npt != 15; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((400+20*$npt))\)
done
for ((npt=0; $npt != 9; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((700+100*$npt))\)
done
