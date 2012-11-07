#!/bin/bash

. setup.bash

hostname
which jsf
echo "#### Start Bases: `date` "
echo "job_date_time=`date +%y%m%d-%H%M%S`"
echo "program_name_version=physsim-110421"

jsf -b -q  bases.C

optstr=
genlumi=`grep GEN_LUMINOSITY run.defs | cut -d"=" -f2`
if [ "x${genlumi}" != "x" -a  -e all.log ] ; then
  xsect=`grep "Total Cross" all.log | cut -d"=" -f2 | cut -d"+" -f1 | sed -e "s/^ *//" `
  if [ "x${xsect}" != "x" -a "x${xsect}" != "xnan" ] ; then
    nevent=`echo "print int(${xsect}*${genlumi})+1" | python -`
    minevents=`grep MAXEVENTS run.defs | cut -d"=" -f2`
    echo "Requested int. luminosity is ${genlumi} fb-1"
    echo "Total cross section of this process is ${xsect} fb"
    echo "Events to be generated from cross section is ${nevent}"
    if [ ${nevent} -lt ${minevents} ] ; then
      nevent=${minevents}
    fi
    optstr=" --maxevt=${nevent} "
  fi
fi
echo "Option for jsf $optstr "

echo "#### Start Spring: `date` "
jsf -b -q --maxevt=${nevent} gui.C

echo "#### End of Spring: `date` "

echo "Create runinfo file."
./mkRuninfo.sh


