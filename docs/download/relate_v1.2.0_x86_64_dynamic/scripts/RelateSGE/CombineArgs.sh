#!/bin/bash

#$ -V
#$ -j y

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

if [[ "$PATH_TO_RELATE" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi

if [ -z ${Ne-} ]
then

  ${PATH_TO_RELATE}/bin/Relate \
    --mode "CombineSections" \
    --chunk_index ${chunk_index} \
    -o ${output} 2>> log/combine_args.log
else
  ${PATH_TO_RELATE}/bin/Relate \
    --mode "CombineSections" \
    -N $Ne \
    --chunk_index ${chunk_index} \
    -o ${output} 2>> log/combine_args.log
fi

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

