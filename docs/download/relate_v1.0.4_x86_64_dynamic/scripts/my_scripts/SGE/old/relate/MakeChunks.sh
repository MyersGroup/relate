#!/bin/bash

#$ -V
#$ -j y
#$ -N convert_from_gp
#$ -P myers.prjc -q short.qc

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

#### check if data exists
check_file_existence (){
  if [ ! -f $1 ]; then
      echo "File " + $1 + " not found!"
      exit 1
  fi
}
check_file_existence sequences.txt
check_file_existence pos.txt
check_file_existence recombination_rate.txt

../bin/./Relate \
  --mode "MakeChunks" \
  --seq sequences.txt \
  --pos pos.txt \
  --rec recombination_rate.txt \
  -n ${chunk_size} 2>> ../log/make_chunks.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0
