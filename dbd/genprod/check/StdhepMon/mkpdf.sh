#!/bin/bash

SEQBGN=106301
SEQEND=106352
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

(
cd log
for i in `seq $SEQBGN $SEQEND` ; do
  pid=${PREF}${i}
  if [ -e ${pid}-part.png ] ; then
    mkpdf ${pid}
    mv -v ${pid}.pdf ../pdf
  fi
done
)

