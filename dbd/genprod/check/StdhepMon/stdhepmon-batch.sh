#!/bin/bash

RUNINFODIR=/home/ilc/miyamoto/DBS/runinfo
PROCIDPREFIX=p
DATADIR=/home/ilc/miyamoto/soft/physsim2/genprod/data
DATADIR=/nfs/g/ilc/soft/samples/gen/CDS/1000
ISHPSS=0
DOANAL=1
MKFIGS=1

# =============== Process just one file.
runone()
{
pid=$1

if [ "x$1" == "x" ] ; then
  echo "Command usage :"
  echo "  $0 [pid] " 
  exit
fi


id
hostname

if [ $ISHPSS -eq 1 ] ; then 
  (
  cd scratch
  getstdhep $pid
  )
  file0=`ls scratch | grep $pid`
  file=scratch/${file0}
else
  file=${DATADIR}/${pid}_01.stdhep
fi

echo $file

if [ $DOANAL -eq 1 ] ; then
  mkdir -p root
  outfile=root/${pid}.root
  jsf -b -q --maxevt=1000000000 --stdhep=${file} \
    --runinfo=${RUNINFODIR} \
    --OutputFile=${outfile} gui.C
  if [ $ISHPSS -eq 1 ] ; then
    rm -vf $file
  fi
fi

if [ $MKFIGS -eq 1 ] ; then
  runinfo=${RUNINFODIR}/${pid}.txt
  procname=`grep process_name ${runinfo} | cut -d"=" -f2 | sed -e "s/ *//g" `
  jsf -b -q "plotfig.C(\"$pid\",\"$procname\")" 
#    rm -f ${outfile}
fi

} 

# ============= main part of the script =====================
# . setup-jsf.sh
. setup.bash
.  /etc/profile.d/hpss.sh
for i in `seq $1 $2 ` ; do 
  pid=${PROCIDPREFIX}${i}
  runinfo=${RUNINFODIR}/${pid}.txt
  if [ -e ${runinfo} ] ; then
    runone ${pid}
  else
    echo "$runinfo does not exist"
  fi
done

