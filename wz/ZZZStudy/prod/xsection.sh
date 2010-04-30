#!/bin/sh
for ((npt=0; $npt != 3; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((270+10*$npt))\)
done
for ((npt=0; $npt != 8; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((300+50*$npt))\)
done
for ((npt=0; $npt != 9; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((700+100*$npt))\)
done
