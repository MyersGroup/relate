#!/bin/bash

#$ -cwd -V
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

chr=${SGE_TASK_ID}
../bin/RelateCoalescentRate --mode ReEstimateBranchLengths -N 30000 -m 1.25e-8 --pos pos_chr${chr}_${pop}.txt --coal relate_${pop}.coal -i relate_chr${chr}_${pop} -o relate_chr${chr}_${pop} 

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

