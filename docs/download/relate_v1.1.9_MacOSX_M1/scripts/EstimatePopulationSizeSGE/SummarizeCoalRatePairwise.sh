#!/bin/bash

#$ -V
#$ -j y
#$ -l hostname=compE*

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

${PATH_TO_RELATE}/bin/RelateCoalescentRate \
  --mode SummarizeCoalescentRateForGenome \
  --chr ${filename_chr} \
  -o ${output}.pairwise

${PATH_TO_RELATE}/bin/RelateCoalescentRate \
  --mode FinalizePopulationSize \
  -o ${output}.pairwise

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

