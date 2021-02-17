#!/bin/bash

#$ -V
#$ -j y
#$ -N convert_from_gp
#$ -P myers.prjc -q short.qc

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

#### check if data exists
check_file_existence (){
  if [ ! -f $1 ]; then
      echo "File " + $1 + " not found!"
      exit 1
  fi
}

check_file_existence ${haps} 
check_file_existence ${sample}
check_file_existence ${dist}
check_file_existence ${map}
check_file_existence ${annot}

${PATH_TO_RELATE}/bin/Relate \
  --mode "MakeChunks" \
  --haps ${haps} \
  --sample ${sample} \
  --map ${map} \
  --dist ${dist} \
  --annot ${annot} \
  --memory ${memory} 2>> log/make_chunks.log

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0
