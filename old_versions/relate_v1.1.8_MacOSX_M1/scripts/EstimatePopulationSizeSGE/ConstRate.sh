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
if [ -z "${bins-}" ];
then

  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
    --mode GenerateConstCoalFile \
    --years_per_gen ${years_per_gen} \
    -i 30000 \
    -o ${output}

else

  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
    --mode GenerateConstCoalFile \
    --years_per_gen ${years_per_gen} \
    --bins ${bins} \
    -i 30000 \
    -o ${output}

fi

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

