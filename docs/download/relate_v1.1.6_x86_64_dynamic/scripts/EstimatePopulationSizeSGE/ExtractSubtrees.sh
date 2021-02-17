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

ind=$((${SGE_TASK_ID}-1))
chromosomes=($(cat ${filename_chr}))
chr=${chromosomes[$ind]}

pop_of_interest=$(echo ${pop_of_interest} | awk '{ gsub(";", ",") ; print $0 }'CC)

${PATH_TO_RELATE}/bin/RelateExtract \
  --mode SubTreesForSubpopulation \
  --pop_of_interest ${pop_of_interest} \
  --poplabels ${filename_poplabels} \
  --anc ${filename}_chr${chr}.anc \
  --mut ${filename}_chr${chr}.mut \
  -o ${output}_${labels}_chr${chr} 2> ${output}_chr${chr}.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

