#!/bin/bash

#$ -V
#$ -j y
#$ -P myers.prjc -q short.qc

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"
 
mkdir est_files
for i in `seq $first_chr 1 $last_chr`
do
  mv relate_chr${i}_${pop}.arg est_files
  mv relate_chr${i}_${pop}.mut est_files
  gunzip relate_chr${i}_${pop}_pre.arg
  mv relate_chr${i}_${pop}_pre.arg relate_chr${i}_${pop}.arg
  cp relate_chr${i}_${pop}_pre.mut relate_chr${i}_${pop}.mut
done
gzip est_files/*

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

