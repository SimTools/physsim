#!/bin/sh
for ((npt=0; $npt != 5; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((150+10*$npt))\)
done
for ((npt=0; $npt != 6; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((200+50*$npt))\)
done
for ((npt=0; $npt != 11; npt=$npt+1))
do 
  jsf -b -q xsection.C\($((500+100*$npt))\)
done
