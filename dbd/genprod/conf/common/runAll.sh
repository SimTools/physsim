#!/bin/bash

. setup.bash

hostname
which jsf
echo "#### Start Bases: `date` "
echo "job_date_time=`date +%y%m%d-%H%M%S`"
echo "program_name_version=physsim-110421"

jsf -b -q  bases.C

echo "#### Start Spring: `date` "
jsf -b -q gui.C

echo "#### End of Spring: `date` "

echo "Create runinfo file."
./mkRuninfo.sh


