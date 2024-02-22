#!/bin/bash

time_dir='/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/time'
for d in `ls -d $time_dir/nnm_100_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done
