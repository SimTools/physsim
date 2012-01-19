#!/bin/bash


mkdir -p log

# bsub -o log/106001.log -J 106001 ./stdhepmon-batch.sh 106001 106003

for ibgn in 106301 106302 106303 `seq 106311 106352` ; do
  ilast=$[${ibgn}+0]
  a=$ibgn
  b=$ilast
  bsub -o log/sub-$a.log -J $a "( ./stdhepmon-batch.sh $a $b > log/$a.log 2>&1 )"
done

