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

first_chunk=0
last_chunk=$(($(ls ${output}_chr${chr}_tmp_chr*.mut | wc -l)-1))

for chunk in `seq ${first_chunk} 1 ${last_chunk}`
do
  cat ${output}_chr${chr}_chr${chunk}.log  
  rm ${output}_chr${chr}_chr${chunk}.log  
done

${PATH_TO_RELATE}/bin/RelateExtract \
  --mode CombineAncMut \
  -o ${output}_chr${chr}_tmp

mv ${output}_chr${chr}_tmp.anc.gz ${output}_chr${chr}.anc.gz
mv ${output}_chr${chr}_tmp.mut.gz ${output}_chr${chr}.mut.gz

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

