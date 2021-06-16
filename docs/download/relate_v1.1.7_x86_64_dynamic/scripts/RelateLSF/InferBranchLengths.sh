#!/bin/bash

#$ -V
#$ -j y

echo "***********************************************"
echo "LSF job ID: "$LSB_JOBID
echo "LSF task ID: "$LSB_JOBINDEX
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

ind=$((${LSB_JOBINDEX}-1))

start_ind=$(($ind * $batch_windows))
end_ind=$((($ind + 1) * $batch_windows - 1))
  
if [ -z ${sample_ages-} ]
then
	if [ -z ${coal-} ]
	then

		if [ -z ${seed-} ]
		then

			${PATH_TO_RELATE}/bin/Relate \
					--mode "InferBranchLengths" \
					-m $mu \
					-N $Ne \
					--first_section $start_ind \
					--last_section $end_ind \
					--chunk_index ${chunk_index} \
					-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		else

			${PATH_TO_RELATE}/bin/Relate \
					--mode "InferBranchLengths" \
					-m $mu \
					-N $Ne \
					--first_section $start_ind \
					--last_section $end_ind \
					--chunk_index ${chunk_index} \
					--seed ${seed} \
					-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		fi

	else

		check_file_existence ${coal} 

		if [ -z ${seed-} ]
		then

			${PATH_TO_RELATE}/bin/Relate \
					--mode "InferBranchLengths" \
					-m $mu \
					--coal ${coal} \
					--first_section $start_ind \
					--last_section $end_ind \
					--chunk_index ${chunk_index} \
					-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		else

			${PATH_TO_RELATE}/bin/Relate \
					--mode "InferBranchLengths" \
					-m $mu \
					--coal ${coal} \
					--first_section $start_ind \
					--last_section $end_ind \
					--chunk_index ${chunk_index} \
					--seed ${seed} \
					-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		fi

	fi
else
	if [ -z ${coal-} ]
	then

		if [ -z ${seed-} ]
		then

			${PATH_TO_RELATE}/bin/Relate \
				--mode "InferBranchLengths" \
				-m $mu \
				-N $Ne \
				--first_section $start_ind \
				--last_section $end_ind \
				--chunk_index ${chunk_index} \
				--sample_ages ${sample_ages} \
				-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		else

			${PATH_TO_RELATE}/bin/Relate \
				--mode "InferBranchLengths" \
				-m $mu \
				-N $Ne \
				--first_section $start_ind \
				--last_section $end_ind \
				--chunk_index ${chunk_index} \
				--seed ${seed} \
				--sample_ages ${sample_ages} \
				-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		fi

	else

		check_file_existence ${coal} 

		if [ -z ${seed-} ]
		then

			${PATH_TO_RELATE}/bin/Relate \
				--mode "InferBranchLengths" \
				-m $mu \
				--coal ${coal} \
				--first_section $start_ind \
				--last_section $end_ind \
				--chunk_index ${chunk_index} \
				--sample_ages ${sample_ages} \
				-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		else

			${PATH_TO_RELATE}/bin/Relate \
				--mode "InferBranchLengths" \
				-m $mu \
				--coal ${coal} \
				--first_section $start_ind \
				--last_section $end_ind \
				--chunk_index ${chunk_index} \
				--seed ${seed} \
				--sample_ages ${sample_ages} \
				-o ${output} 2>> ${LSB_JOBINDEX}_infer_branch_length_c${chunk_index}.log 

		fi

	fi
fi

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

