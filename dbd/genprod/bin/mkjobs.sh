#!/bin/bash

# =====================================================================
Initialize()
{
  TOPDIR=${KFLIBROOT}/dbd
  JOBDIR=${TOPDIR}/jobs
  _MAXEVENTS=50000

  echo "mkjob.sh : creates directories to run physsim jobs."
  echo "All path should be given in absolute path. Otherwise it may not run correctly."
  read -p "Directory where genprod/ and TTHStudy/ etc can be found ? [${TOPDIR}] : " topdir
  if [ "x$topdir" != "x" ] ; then TOPDIR=${topdir} ; fi
  JOBDIR=`pwd`/jobs
  read -p "Job directory ? [${JOBDIR}] : " jobdir
  if [ "x${jobdir}" != "x" ] ; then JOBDIR=${jobdir} ; fi

  CONFDIR=${TOPDIR}/genprod/conf
  processlist=${CONFDIR}/common/process.list
  read -p "Process list file ? [${processlist}] : " processlist_in
  if [ "x${processlist_in}" != "x" ] ; then processlist=${processlist_in} ; fi
  
  read -p "Number of events to generate ? [${_MAXEVENTS}] : " maxevents
  if [ "x${maxevents}" != "x" ] ; then _MAXEVENTS=${maxevents} ; fi

  echo " ------------------------------------------ "
  echo " Top directory : ${TOPDIR}"
  echo " Job directory : ${JOBDIR}"
  echo " Conf directory : ${CONFDIR}"
  echo " Process list file : ${processlist}"
  echo " Max number of events to generate : ${_MAXEVENTS} "
  read -p "OK to continue ? (CTRL-C to break)" ans
  
}


# =============================================
mkrunfiles()
{
  local sedfile=temp$$.sed
  rm -f ${sedfile}
  for file in $1 $2 $3 ; do
    while read line ; do
      local varname=`echo $line | cut -d'=' -f1`
      local value=`echo $line | cut -d'=' -f2`
      echo "s|@@${varname}@@|${value}|g" >> ${sedfile}
    done < ${file}
  done

  for f in *.in ; do 
    newf=`echo ${f} | sed -e "s/\.in$//" `
    sed -f ${sedfile} ${f} > ${newf}
    rm -f ${f}
  done

  chmod +x *.sh
  if [ -e prepfiles.sh ] ; then 
    . prepfiles.sh 
    rm -f prepfiles.sh
  fi
  rm -f ${sedfile}
}


#====================================================
polstr()
{
cat <<EOF | python - $1
import sys
argvs=sys.argv
polstr="R"
epol=float(argvs[1])
if epol < 0.0 :
  polstr="L"
polval=int(abs(epol)*100.0)
if polval > 99 :
  print polstr
else :
  print polstr+str(polval)
EOF

}

# =============================================
mkrundef()
{
  out=$1
  echo "ProcessID=${pid}" > ${out}
  echo "PolElectron=${epol}" >> ${out}
  echo "PolPositron=${ppol}" >> ${out}

  beampara="E1000-Aug312010"
  genvers="Gp01-02"
  epolstr=e`polstr ${epol}`
  ppolstr=p`polstr ${ppol}`

  case "${ttmodes}" in 
    "all") echo "WmModesLo=1" >> ${out} ;
	 echo "WmModesHi=12" >> ${out} ;
	 echo "WpModesLo=1" >> ${out} ;
         echo "WpModesHi=12" >> ${out} ;
         echo "FinalStateMix=0" >> ${out} ;;
    "6q") echo "WmModesLo=4" >> ${out} ;
	 echo "WmModesHi=12" >> ${out} ;
	 echo "WpModesLo=4" >> ${out} ;
         echo "WpModesHi=12" >> ${out} ;
         echo "FinalStateMix=0" >> ${out} ;;
    "ln4q") echo "WmModesLo=1" >> ${out} ;
	 echo "WmModesHi=3" >> ${out} ;
	 echo "WpModesLo=4" >> ${out} ;
         echo "WpModesHi=12" >> ${out} ;
         echo "FinalStateMix=1" >> ${out} ;;
    "en4q") echo "WmModesLo=1" >> ${out} ;
	 echo "WmModesHi=1" >> ${out} ;
	 echo "WpModesLo=4" >> ${out} ;
         echo "WpModesHi=12" >> ${out} ;
         echo "FinalStateMix=1" >> ${out} ;;
    "2l2nbb") echo "WmModesLo=1" >> ${out} ;
	 echo "WmModesHi=3" >> ${out} ;
	 echo "WpModesLo=1" >> ${out} ;
         echo "WpModesHi=3" >> ${out} ;
         echo "FinalStateMix=0" >> ${out} ;;
    *)
      echo "Fatal error : ttmodes ${ttmodes} is not defined"
      exit -1
  esac

  
  case "${hmodes}" in
    "zqq") echo "ZModesLo=7" >> ${out} ;  echo "ZModesHi=12" >> ${out} ;
	   echo "DecayModesForH=0" >> ${out} ;;
    "zll") echo "ZModesLo=4" >> ${out} ;  echo "ZModesHi=6" >> ${out} ;
	   echo "DecayModesForH=0" >> ${out} ;;
    "znn") echo "ZModesLo=1" >> ${out} ;  echo "ZModesHi=3" >> ${out} ;
	   echo "DecayModesForH=0" >> ${out} ;;
    "gbb") echo "DecayModesForH=0" >> ${out} ;;
    "hbb") echo "DecayModesForH=5" >> ${out} ;;
    "hnonbb") echo "DecayModeForH=105" >> ${out} ;;
    "all") echo "DecayModeForH=0" >> ${out} ;
           echo "ZModesLo=1" >> ${out} ; echo "ZModesHi=12" >> ${out} ;;
    *)
      echo "Fatal error : hmodes ${hmodes} is not defined."
      exit -1
   esac

   prstr=${process}-${ttmodes}-${hmodes}
   echo "StdhepTitle=physsim-${prstr}" >> ${out}
#   echo "StdhepFileName=p${pid}" >> ${out}
   echo "StdhepFileName=${beampara}.P${prstr}.${epolstr}.${ppolstr}.${genvers}.I${pid}"
   echo "MAXEVENTS=${_MAXEVENTS}" >> ${out}
}

# =============================================
mkfiles()
{
  pid=$1
  Process=$2
  Ttmodes=$3
  Hmodes=$4
  epol=$5
  ppol=$6

  process=`echo $Process | tr [:upper:] [:lower:]`
  hmodes=`echo $Hmodes | tr [:upper:] [:lower:]`
  ttmodes=`echo $Ttmodes | tr [:upper:] [:lower:]`

  processstr=p${pid}
  rundir=${JOBDIR}/${processstr}
  origdir=`pwd`

  if [ -e ${rundir} ] ; then
#    rm -rf ${rundir}
    echo "Output directory, ${rundir}, exist already."
    echo "Delete output directory first."
    exit
  fi

  echo "Creating run directory, ${rundir}"
  mkdir -p ${rundir}
  cd ${rundir}

  cp -rp ${CONFDIR}/common/* ./  
  cat ${CONFDIR}/${process}/process.conf.in \
      ${CONFDIR}/common/common.conf.in > jsf.conf.in
  cp -p  ${CONFDIR}/dbs/common.defs  runinfo_common.defs
  mkrundef run.defs 
  mkrunfiles ${CONFDIR}/${process}/process.defs ${CONFDIR}/common/common.defs run.defs

  cd ${origdir}

}


Initialize

while read line ; do
  if [ "x${line:0:1}" == "x#" -o "x${line:0:1}" == "x" ] ; then continue ; fi
  echo "mkfiles ${line}"
  mkfiles ${line}
done < ${processlist}

( 
cd ${JOBDIR} 
  echo "for d in p* ; do ( cd \${d} && . subAll.sh ) ; done " > allsub.sh 
)

