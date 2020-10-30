#!/bin/bash

#$ -V
#$ -j y

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

if [[ "$PATH_TO_RELATE" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi

${PATH_TO_RELATE}/bin/Relate \
  --mode "FindEquivalentBranches" \
  --chunk_index ${chunk_index} \
  -o ${output} 2>> log/find_equivalent_branches_c${chunk_index}.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

