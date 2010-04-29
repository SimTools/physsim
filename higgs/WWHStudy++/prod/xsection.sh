#!/bin/sh
for ((npt=0; $npt != 5; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((1100+100*$npt))\)
done
exit 0
for ((npt=0; $npt != 6; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((280+20*$npt))\)
done
for ((npt=0; $npt != 4; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((400+50*$npt))\)
done
for ((npt=0; $npt != 5; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((600+100*$npt))\)
done
