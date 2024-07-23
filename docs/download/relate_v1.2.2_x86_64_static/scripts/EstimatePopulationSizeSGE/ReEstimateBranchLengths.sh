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

chunk=$((${SGE_TASK_ID}-1))

if [ -z "${seed-}" ];
then

  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
    --mode ReEstimateBranchLengths \
    --coal ${output}.coal \
    --dist ${output}_chr${chr}.dist \
    -m ${mu} \
    -i ${output}_chr${chr}_tmp_chr${chunk} \
    -o ${output}_chr${chr}_tmp_chr${chunk} 2> ${output}_chr${chr}_chr${chunk}.log 

else

  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
    --mode ReEstimateBranchLengths \
    --coal ${output}.coal \
    --dist ${output}_chr${chr}.dist \
    --seed $seed \
    -m ${mu} \
    -i ${output}_chr${chr}_tmp_chr${chunk} \
    -o ${output}_chr${chr}_tmp_chr${chunk} 2> ${output}_chr${chr}_chr${chunk}.log 
fi

rm ${output}_chr${chr}_tmp_chr${chunk}.anc.gz
rm ${output}_chr${chr}_tmp_chr${chunk}.mut.gz

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

