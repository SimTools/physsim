#!/bin/bash

SEQBGN=106401
SEQEND=106452
PREF=

mkpdf()
{
  pid=$1
  flist=""
  for t in title part event1 event2 2jets 4jets ; do 
    convert ${pid}-${t}.png ${pid}-${t}.pdf 
    flist="${flist} ${pid}-${t}.pdf"
  done
  pdftk ${flist} output ${pid}.pdf
  rm -f ${flist}
}

mkdir -p pdf
filelist=`pwd`/filelist.txt

(
cd log
for i in `seq $SEQBGN $SEQEND` ; do
  pid=${PREF}${i}
  echo $pid
  if [ -e ${pid}-part.png ] ; then
    mkpdf ${pid}
    newname=`grep I${pid} $filelist`
    mv -v ${pid}.pdf ../pdf/${newname}.pdf
  fi
done
)

