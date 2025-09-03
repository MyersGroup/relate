#!/bin/bash

#$ -V
#$ -j y

echo "***********************************************"
echo "LSF job ID: "$LSB_JOBID
echo "LSF task ID: "$LSB_JOBINDEX
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

${PATH_TO_RELATE}/bin/Relate \
  --mode "FindEquivalentBranches" \
  --chunk_index ${chunk_index} \
  -o ${output} 2>> log/find_equivalent_branches_c${chunk_index}.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

