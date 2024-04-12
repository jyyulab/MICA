#!/bin/bash

dir='/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/computing_performance/pecanpy/hyperparameter'
for d in `ls -d $dir/phyper_1.8_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done

for d in `ls -d $dir/phyper_2.0_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done

for d in `ls -d $dir/phyper_2.2_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done

for d in `ls -d $dir/phyper_2.4_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done

for d in `ls -d $dir/phyper_2.6_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done

for d in `ls -d $dir/phyper_2.8_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done

for d in `ls -d $dir/phyper_3.0_*`
do
  echo $d/MICA_GE.sh
  bsub < $d/MICA_GE.sh
done
