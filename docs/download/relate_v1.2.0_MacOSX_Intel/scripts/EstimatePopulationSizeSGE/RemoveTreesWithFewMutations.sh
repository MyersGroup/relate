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

#delete all trees that have fewer than $threshold mutations
${PATH_TO_RELATE}/bin/RelateExtract \
  --mode RemoveTreesWithFewMutations \
  --threshold $threshold \
  --anc ${filename}_chr${chr}.anc \
  --mut ${filename}_chr${chr}.mut \
  -o ${output}_chr${chr} 2> ${output}_chr${chr}.log

gzip ${output}_chr${chr}.anc
gzip ${output}_chr${chr}.mut

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

