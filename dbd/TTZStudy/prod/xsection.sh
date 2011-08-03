#!/bin/sh
for ((npt=0; $npt != 57; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((440+10*$npt))\)
done
