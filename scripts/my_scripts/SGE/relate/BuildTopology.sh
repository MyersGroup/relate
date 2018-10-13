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

ind=$((${SGE_TASK_ID}-1))

start_ind=$(($ind * $batch_windows))
end_ind=$((($ind + 1) * $batch_windows - 1))
  
../bin/./Relate \
    --mode "BuildTopology" \
    --first_section $start_ind \
    --last_section $end_ind \
    --chunk_index ${chunk_index} \
    -o ${output} 2>> ${SGE_TASK_ID}_build_c${chunk_index}.log 

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

