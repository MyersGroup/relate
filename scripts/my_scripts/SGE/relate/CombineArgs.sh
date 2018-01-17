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

../bin/./Relate \
  --mode "CombineSections" \
  -N $Ne \
  --chunk_index ${chunk_index} \
  -o ${output} 2>> log/combine_args.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

