#!/bin/sh
for ((npt=0; $npt != 16; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((250+50*$npt))\)
done
