#!/bin/sh
for ((npt=0; $npt != 54; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((470+10*$npt))\)
done
