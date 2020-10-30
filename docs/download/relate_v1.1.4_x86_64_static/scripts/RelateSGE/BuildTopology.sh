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

ind=$((${SGE_TASK_ID}-1))

start_ind=$(($ind * $batch_windows))
end_ind=$((($ind + 1) * $batch_windows - 1))

if [[ "$PATH_TO_RELATE" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi

painting=$(echo ${painting} | awk '{ gsub(";", ",") ; print $0 }'CC)
if [ -z ${seed-} ]			 
then
	${PATH_TO_RELATE}/bin/Relate \
			--mode "BuildTopology" \
			--first_section $start_ind \
			--last_section $end_ind \
			--chunk_index ${chunk_index} \
			--painting ${painting} \
			-o ${output} 2>> ${SGE_TASK_ID}_build_c${chunk_index}.log 
else
	${PATH_TO_RELATE}/bin/Relate \
		--mode "BuildTopology" \
		--first_section $start_ind \
		--last_section $end_ind \
		--chunk_index ${chunk_index} \
		--seed ${seed} \
		--painting ${painting} \
		-o ${output} 2>> ${SGE_TASK_ID}_build_c${chunk_index}.log 
fi
echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

