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

check_file_existence ${haps} 
check_file_existence ${sample}
check_file_existence ${map}

if [[ "$PATH_TO_RELATE" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi


if [ ! -z ${dist-} ]
then
  check_file_existence ${dist}

  if [ -z ${annot-} ]
  then

    ${PATH_TO_RELATE}/bin/Relate \
      --mode "MakeChunks" \
      --haps ${haps} \
      --sample ${sample} \
      --map ${map} \
      --dist ${dist} \
      --memory ${memory} \
      --output $output 2>> log/make_chunks.log

  else

    check_file_existence ${annot}
    ${PATH_TO_RELATE}/bin/Relate \
      --mode "MakeChunks" \
      --haps ${haps} \
      --sample ${sample} \
      --map ${map} \
      --dist ${dist} \
      --annot ${annot} \
      --memory ${memory} \
      --output $output 2>> log/make_chunks.log

  fi
else
  if [ -z ${annot-} ]
  then

    ${PATH_TO_RELATE}/bin/Relate \
      --mode "MakeChunks" \
      --haps ${haps} \
      --sample ${sample} \
      --map ${map} \
      --memory ${memory} \
      --output $output 2>> log/make_chunks.log

  else

    check_file_existence ${annot}
    ${PATH_TO_RELATE}/bin/Relate \
      --mode "MakeChunks" \
      --haps ${haps} \
      --sample ${sample} \
      --map ${map} \
      --annot ${annot} \
      --memory ${memory} \
      --output $output 2>> log/make_chunks.log

  fi
fi



echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0
