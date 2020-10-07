#!/bin/bash

#$ -V
#$ -j y

echo "***********************************************"
echo "Slurm job ID: "$SLURM_JOBID
echo "Slurm task ID: "$SLURM_ARRAY_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

if [[ "$PATH_TO_RELATE" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi
painting=$(echo ${painting} | awk '{ gsub(";", ",") ; print $0 }'CC)
## paint all sequences against each other
${PATH_TO_RELATE}/bin/Relate --mode "Paint" --chunk_index ${chunk_index} --painting ${painting} --output ${output} 2>> log/paint_c${chunk_index}.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

