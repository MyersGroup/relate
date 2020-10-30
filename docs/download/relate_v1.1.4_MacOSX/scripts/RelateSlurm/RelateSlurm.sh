#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/RelateSlurm/RelateSlurm.sh"
  echo ""
  echo "--haps:   Filename of haps file (Output file format of Shapeit)."
  echo "--sample: Filename of sample file (Output file format of Shapeit)."
  echo "--map:    Genetic map."
  echo "-m,--mu:  Mutation rate."
  echo "-N,--Ne:  Effective population size."
  echo "-o,--output: Filename of output without file extension."
  echo "--slurm_options: All slurm-specific arguments, e.g., account."
  echo "--dist:   Optional but recommended. Distance in BP between SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used."
  echo "--annot:  Optional. Filename of file containing additional annotation of snps. Can be generated using RelateFileFormats."
  echo "--memory: Optional. Approximate memory allowance in GB passed to Relate for storing distance matrices. Default is 5GB."
	echo "--sample_ages: Optional. Filename of file containing sample ages (one per line)." 
	echo "--coal:   Optional. Filename of file containing coalescent rates. If specified, it will overwrite --effectiveN."
  echo "--painting: Optional. Copying and transition parameters in
                    chromosome painting algorithm. Format: theta,rho. Default: 0.025,1." 
  echo "--seed:   Optional. Seed for MCMC in branch lengths estimation."
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/RelateSlurm/RelateSlurm.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

memory=5
painting="0.025,1"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --haps)
      haps="$2"
      shift # past argument
      shift # past value
      ;;
    --sample)
      sample="$2"
      shift # past argument
      shift # past value
      ;;
    --map)
      map="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--mu)
      mu="$2"
      shift # past argument
      shift # past value
      ;;
    -N|--Ne)
      Ne="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      output="$2"
      shift # past argument
      shift # past value
      ;;
    --seed)
      seed="$2"
      shift # past argument
      shift # past value
      ;;  
    --dist)
      dist="$2"
      shift # past argument
      shift # past value
      ;;
    --annot)
      annot="$2"
      shift # past argument
      shift # past value
      ;;
    --memory)
      memory="$2"
      shift # past argument
      shift # past value
      ;;
		--sample_ages)
			sample_ages="$2"
			shift # past argument
			shift # past value
			;;
		--painting)
			painting="$2"
			shift # past argument
			shift # past value
			;;
    --coal)
      coal="$2"
      shift # past argument
      shift # past value
      ;;
    --slurm_options)
      slurm_options="$2"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

if [ -z ${mem-} ];
then
  mem=$((2*${memory}))
fi

echo "********************************"
echo "Parameters passed to Relate:"
echo "haps   = $haps"
echo "sample = $sample"
echo "map    = $map"
if [ ! -z "${coal-}" ];
then
  echo "coal   = $coal"
else 
  echo "Ne     = $Ne"
fi
echo "mu     = $mu"
echo "output = $output"
if [ ! -z "${dist-}" ];
then
  echo "dist   = $dist"
fi
if [ ! -z "${annot-}" ];
then
  echo "annot  = $annot"
fi
if [ ! -z "${memory-}" ];
then
  echo "memory = $memory"
fi
if [ ! -z "${sample_ages-}" ];
then
	echo "sample_ages = $sample_ages"
fi
if [ ! -z "${seed-}" ];
then
  echo "seed   = $seed"
fi
if [ ! -z "${slurm_options-}" ];
then
  echo "slurm_options = $slurm_options"
fi
echo "********************************"

######################################################################################################

#### create dir for logfiles
if [ -d "${output}" ]; then
  # Control will enter here if $DIRECTORY exists.
  echo ""
  echo "Warning: Directory ${output} exists. Strongly recommended to run this script with a unique output filename or after removing this directory."
  echo ""
fi
mkdir -p ${output}
mkdir -p ${output}/log
cd ${output}

######################################################################################################

#### DO NOT CHANGE THE FOLLOWING PARAMETERS
batch_windows=5

######################################################################################################

#### Divide data into chunks of 500k SNPs (the chunks are overlapping to avoid edge effects)
#-sync y causes sbatch to wait for the job to complete before exiting.

if [[ "$haps" != /* ]]; 
then
  haps="../${haps}"
fi
if [[ "$sample" != /* ]]; 
then
  sample="../${sample}"
fi
if [[ "$map" != /* ]]; 
then
  map="../${map}"
fi

if [ ! -z ${dist-} ]
then

  if [[ "$dist" != /* ]]; 
  then
    dist="../${dist}"
  fi

  if [ -z ${annot-} ]
  then

    jid=$(sbatch -J "make_chunks_${output}" \
         -W \
         --parsable \
         --error=log/make_chunks.log \
         --output=log/make_chunks.log \
         --export PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},dist=${dist},output=${output} \
         ${slurm_options} \
         "${PATH_TO_RELATE}/scripts/RelateSlurm/MakeChunks.sh")

  else

    if [[ "$annot" != /* ]]; 
    then
      annot="../${annot}"
    fi

    jid=$(sbatch -J "make_chunks_${output}" \
         -W \
         --parsable \
         --error=log/make_chunks.log \
         --output=log/make_chunks.log \
         --export PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},dist=${dist},annot=${annot},output=${output} \
         ${slurm_options} \
         "${PATH_TO_RELATE}/scripts/RelateSlurm/MakeChunks.sh")

  fi
else
  if [ -z ${annot-} ]
  then

    jid=$(sbatch -J "make_chunks_${output}" \
         -W \
         --parsable \
         --error=log/make_chunks.log \
         --output=log/make_chunks.log \
         --export=PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},output=${output} \
         ${slurm_options} \
         "${PATH_TO_RELATE}/scripts/RelateSlurm/MakeChunks.sh")
 
  else

    if [[ "$annot" != /* ]]; 
    then
      annot="../${annot}"
    fi

    jid=$(sbatch -J "make_chunks_${output}" \
         -W \
         --parsable \
         --error=log/make_chunks.log \
         --output=log/make_chunks.log \
         --export PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},annot=${annot},output=${output} \
         ${slurm_options} \
         "${PATH_TO_RELATE}/scripts/RelateSlurm/MakeChunks.sh")

  fi
fi


###########GET NUMBER OF CHUNKS
num_chunks=($(${PATH_TO_RELATE}/scripts/RelateSlurm/read_bin.py "${output}/parameters.bin"))
num_chunks=${num_chunks[2]}
prev_chunk=-1

echo "********************************"
echo "Number of chunks: "$num_chunks

######################################################################################################
#### Build trees for each chunk
painting=$(echo ${painting} | awk '{ gsub(",", ";") ; print $0 }'CC)
prev_jid=(${jid})
num_paintings=3

## paint all sequences against each other
for chunk in `seq 0 $(($num_chunks - 1))`;
do

  #GET NUMBER OF WINDOWS
  num_windows=($(${PATH_TO_RELATE}/scripts/RelateSlurm/read_bin.py "${output}/parameters_c"${chunk}".bin"))
  num_windows=${num_windows[2]}
  num_batched_windows=$(($num_windows/$batch_windows + 1))

  wait_jid=${prev_jid[0]}
  if [ ${#prev_jid[@]} -gt ${num_paintings} ];
  then
    wait_jid=${prev_jid[$(($chunk - ${num_paintings} + 1))]}
  fi

  ## paint all sequences against each other
  #make sure that only 5 paintings exist at a time. Chunk 5 is painted only after chunk 0 is done etc.
  jid=$(sbatch --depend afterok:${wait_jid} \
       --parsable \
       -J paint_${output}_${chunk} \
       --export PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=${chunk},painting=${painting},output=${output} \
       -e log/paint_c${chunk}.log \
       -o log/paint_c${chunk}.log \
       ${slurm_options} \
       ${PATH_TO_RELATE}/scripts/RelateSlurm/Paint.sh)
				 
	if [ -z ${seed-} ]			 
	then
		## build tree topologies
		jid=$(sbatch --depend afterok:${jid} \
				 --parsable \
				 -J build_topology_${output}_${chunk} \
				 --array 1-$num_batched_windows \
				 --export PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output},batch_windows=$batch_windows,painting=${painting} \
				 -e build_${output}.log \
				 -o build_${output}.log \
				 ${slurm_options} \
				 ${PATH_TO_RELATE}/scripts/RelateSlurm/BuildTopology.sh)
	else
		## build tree topologies
		jid=$(sbatch --depend afterok:${jid} \
			--parsable \
			-J build_topology_${output}_${chunk} \
			--array 1-$num_batched_windows \
			--export PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output},batch_windows=$batch_windows,seed=${seed},painting=${painting} \
			-e build_${output}.log \
			-o build_${output}.log \
			${slurm_options} \
			${PATH_TO_RELATE}/scripts/RelateSlurm/BuildTopology.sh)
	fi

  ## find equivalent branches in adjacent trees 
  jid=$(sbatch --depend afterok:${jid} \
       --parsable \
       -J find_equivalent_branches_${output}_${chunk} \
       --export PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output} \
       -e log/find_equivalent_branches_c${chunk}.log \
       -o log/find_equivalent_branches_c${chunk}.log \
       ${slurm_options} \
       ${PATH_TO_RELATE}/scripts/RelateSlurm/FindEquivalentBranches.sh)

  prev_jid+=(${jid})

  ## infer branch lengths
  if [ -z ${sample_ages-} ]
	then
		if [ -z ${coal-} ]
		then
			if [ -z ${seed-} ]
			then
				jid=$(sbatch --depend afterok:${jid} \
             --parsable \
						 -J infer_branch_lengths_${output}_${chunk} \
						 --array 1-$num_batched_windows \
						 --export PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
             ${slurm_options} \
             ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			else
				jid=$(sbatch --depend afterok:${jid} \
             --parsable \
						 -J infer_branch_lengths_${output}_${chunk} \
						 --array 1-$num_batched_windows \
						 --export PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output},seed=${seed} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
             ${slurm_options} \
             ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			fi

			## combine args into one file
			jid=$(sbatch --depend afterok:${jid} \
           --parsable \
					 -J combine_args_${output} \
					 --export PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,chunk_index=${chunk},output=${output} \
					 -e log/combine_args_c${chunk}.log \
					 -o log/combine_args_c${chunk}.log \
           ${slurm_options} \
           ${PATH_TO_RELATE}/scripts/RelateSlurm/CombineArgs.sh)

		else

			if [[ "$coal" != /* ]]; 
			then
				coal="../${coal}"
			fi

			if [ -z ${seed-} ]
			then
				jid=$(sbatch --depend afterok:${jid} \
             --parsable \
						 -J infer_branch_lengths_${output}_${chunk} \
						 --array 1-$num_batched_windows \
						 --export PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
             ${slurm_options} \
             ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			else
				jid=$(sbatch --depend afterok:${jid} \
             --parsable \
						 -J infer_branch_lengths_${output}_${chunk} \
						 --array 1-$num_batched_windows \
						 --export PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output},seed=${seed} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
             ${slurm_options} \
             ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			fi

			## combine args into one file
			jid=$(sbatch --depend afterok:${jid} \
           --parsable \
					 -J combine_args_${output} \
					 --export PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},chunk_index=${chunk},output=${output} \
					 -e log/combine_args_c${chunk}.log \
					 -o log/combine_args_c${chunk}.log \
           ${slurm_options} \
           ${PATH_TO_RELATE}/scripts/RelateSlurm/CombineArgs.sh)

		fi

  else

		if [ -z ${coal-} ]
		then
			if [ -z ${seed-} ]
			then
				jid=$(sbatch --depend afterok:${jid} \
          --parsable \
					-J infer_branch_lengths_${output}_${chunk} \
					--array 1-$num_batched_windows \
					--export PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
          ${slurm_options} \
          ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			else
				jid=$(sbatch --depend afterok:${jid} \
          --parsable \
					-J infer_branch_lengths_${output}_${chunk} \
					--array 1-$num_batched_windows \
					--export PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output},seed=${seed} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
          ${slurm_options} \
          ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			fi

			## combine args into one file
			jid=$(sbatch --depend afterok:${jid} \
        --parsable \
				-J combine_args_${output} \
				--array ${PWD}/${output} \
				--export PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,chunk_index=${chunk},output=${output} \
				-e log/combine_args_c${chunk}.log \
				-o log/combine_args_c${chunk}.log \
        ${slurm_options} \
        ${PATH_TO_RELATE}/scripts/RelateSlurm/CombineArgs.sh)

		else

			if [[ "$coal" != /* ]]; 
			then
				coal="../${coal}"
			fi

			if [ -z ${seed-} ]
			then
				jid=$(sbatch --depend afterok:${jid} \
          --parsable \
					-J infer_branch_lengths_${output}_${chunk} \
					--array 1-$num_batched_windows \
					--export PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
          ${slurm_options} \
          ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			else
				jid=$(sbatch --depend afterok:${jid} \
          --parsable \
					-J infer_branch_lengths_${output}_${chunk} \
					--array 1-$num_batched_windows \
					--export PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output},seed=${seed} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
          ${slurm_options} \
          ${PATH_TO_RELATE}/scripts/RelateSlurm/InferBranchLengths.sh)
			fi

			## combine args into one file
			jid=$(sbatch --depend afterok:${jid} \
        --parsable \
				-J combine_args_${output} \
				--export PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},chunk_index=${chunk},output=${output} \
				-e log/combine_args_c${chunk}.log \
				-o log/combine_args_c${chunk}.log \
        ${slurm_options} \
        ${PATH_TO_RELATE}/scripts/RelateSlurm/CombineArgs.sh)

		fi
	fi

  prev_chunk=$chunk

done

############################# Finalize ###########################
#-sync y causes sbatch to wait for the job to complete before exiting. 
#finalize results
if [ -z ${sample_ages-} ]
then
	jid=$(sbatch -W \
       --parsable \
			 --depend afterok:${jid} \
			 -J finalize_${output} \
			 --export PATH_TO_RELATE=${PATH_TO_RELATE},output=${output} \
			 -e log/combine_args.log \
			 -o log/combine_args.log \
       ${slurm_options} \
       ${PATH_TO_RELATE}/scripts/RelateSlurm/Finalize.sh)
else
	jid=$(sbatch -W \
    --parsable \
		--depend afterok:${jid} \
		-J finalize_${output} \
		--export PATH_TO_RELATE=${PATH_TO_RELATE}, sample_ages=${sample_ages},output=${output} \
		-e log/combine_args.log \
		-o log/combine_args.log \
    ${slurm_options} \
    ${PATH_TO_RELATE}/scripts/RelateSlurm/Finalize.sh)
fi
##################################################################
#clean up directory
mv *.log log/
tar -cf log.tar log/
rm -rf log/

gzip ${output}.anc
gzip ${output}.mut

cd ..

