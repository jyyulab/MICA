#!/bin/bash

dir='/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/computing_performance/EU-kNN'
for d in `ls -d $dir/nne_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done
