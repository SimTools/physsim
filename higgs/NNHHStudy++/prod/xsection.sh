#!/bin/sh
for ((npt=0; $npt != 23; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((400+50*$npt))\)
done
