#!/bin/bash

echo "***********************************************"
echo "Slurm job ID: "$SLURM_JOBID
echo "Slurm task ID: "$SLURM_ARRAY_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

ind=$((${SLURM_ARRAY_TASK_ID}-1))

start_ind=$(($ind * $batch_windows))
end_ind=$((($ind + 1) * $batch_windows - 1))

if [[ "$PATH_TO_RELATE" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi

painting=$(echo ${painting} | awk '{ gsub(";", ",") ; print $0 }'CC)

args=()

args+=( "--mode BuildTopology" )
args+=( "--chunk_index ${chunk_index}" )
args+=( "--first_section $start_ind" )
args+=( "--last_section $end_ind" )
args+=( "--painting ${painting}" )
args+=( "-o ${output}" )
[[ ! -z "${sample_ages-}" ]] && args+=( "--sample_ages ${sample_ages}" )
[[ ! -z "${seed-}" ]] && args+=( "--seed ${seed}" )
[[ ${fb-} -ne -1 ]] && args+=( "--fb ${fb}" )
[[ ! -z "${consistency-}" ]] && args+=( "--no_consistency" )

${PATH_TO_RELATE}/bin/Relate ${args[@]} 2>> ${SLURM_ARRAY_TASK_ID}_build_c${chunk_index}.log 


echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

