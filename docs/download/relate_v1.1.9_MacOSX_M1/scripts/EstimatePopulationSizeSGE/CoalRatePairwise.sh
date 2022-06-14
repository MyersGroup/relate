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

bins=$(echo ${bins} | awk '{ gsub(";", ",") ; print $0 }'CC)
ind=$((${SGE_TASK_ID}-1))
chromosomes=($(cat ${filename_chr}))
chr=${chromosomes[$ind]}

if [ -z "${bins-}" ];
then

  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
    --mode CoalescentRateForSection \
    --years_per_gen ${years_per_gen} \
    --dist ${output}_chr${chr}.dist \
    --chr ${filename_chr} \
    -i ${output}_chr${chr} \
    -o ${output}.pairwise_chr${chr}

else

  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
    --mode CoalescentRateForSection \
    --years_per_gen ${years_per_gen} \
    --bins ${bins} \
    --dist ${output}_chr${chr}.dist \
    -i ${output}_chr${chr} \
    -o ${output}.pairwise_chr${chr}

fi

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

