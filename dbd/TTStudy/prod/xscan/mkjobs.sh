#!/bin/bash


min_droots=-10
topdir=`pwd`
jobtop=$topdir/jobs-QCDISRBS
scripts=${topdir}/scripts
subsh=sub.sh

if [ -e ${jobtop}/suball.sh ] ; then 
  mv ${jobtop}/suball.sh ${jobtop}/suball.sh-`datestr`
fi

for i in `seq 0 40` ; do 
# for i in `seq 0 100` ; do 
  jobserstr=`printf %3.3d $i`
  jobdir=${jobtop}/job${jobserstr}
  mkdir -p ${jobdir}
  (
    cd ${jobdir}
    droots=`echo print "${min_droots}+0.5*$i" | python -`
    sed -e "s/%%DELTA_ROOTS%%/${droots}/" ${scripts}/jsf.conf.in > jsf.conf
    ln -s ${scripts}/TTSpring.so .
    ln -s ${scripts}/bases.C .
    ln -s ${scripts}/setup.bash .
    echo "#!/bin/bash" > ${subsh}
    echo "bsub -o bases.log -q e -J job${jobserstr} \"( . setup.bash && jsf -b -q bases.C > bases.log 2>&1 )\" " > ${subsh}
  )
  echo " ( cd job${jobserstr} && . ${subsh} )" >> ${jobtop}/suball.sh 
done
echo "Jobs have been written in ${jobtop}"

