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

## paint all sequences against each other
painting=$(echo ${painting} | awk '{ gsub(";", ",") ; print $0 }'CC)
${PATH_TO_RELATE}/bin/Relate --mode "Paint" --chunk_index ${chunk_index} --output ${output} --painting ${painting} 2>> log/paint_c${chunk_index}.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

