#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/RelateLSF/RelateLSF.sh"
  echo ""
  echo "--haps:   Filename of haps file (Output file format of Shapeit)."
  echo "--sample: Filename of sample file (Output file format of Shapeit)."
  echo "--map:    Genetic map."
  echo "-m,--mu:  Mutation rate."
  echo "-N,--Ne:  Effective population size."
  echo "-o,--output: Filename of output without file extension."
  echo "-P:       Assign job to the project."
  echo "-q:       Queue."
  echo "--min_mem: Optional. Min_mem."
  echo "--max_mem: Optional. Max_mem."
  echo "--dist:   Optional but recommended. Distance in BP between SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used."
  echo "--annot:  Optional. Filename of file containing additional annotation of snps. Can be generated using RelateFileFormats."
  echo "--memory: Optional. Approximate memory allowance in GB for storing distance matrices. Default is 5GB."
	echo "--sample_ages: Optional. Filename of file containing sample ages (one per line)." 
	echo "--coal:   Optional. Filename of file containing coalescent rates. If specified, it will overwrite --effectiveN." 
	echo "--painting: Optional. Copying and transition parameters in
                  	chromosome painting algorithm. Format: theta,rho. Default: 0.025,1." 
  echo "--seed:   Optional. Seed for MCMC in branch lengths estimation."
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/RelateLSF/RelateLSF.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

min_mem="5000"
max_mem="20000"
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
    --coal)
      coal="$2"
      shift # past argument
      shift # past value
      ;;
		--painting)
			painting="$2"
			shift # past argument
			shift # past value
			;;
    -P)
      p="$2"
      shift # past argument
      shift # past value
      ;;
    -q)
      q="$2"
      shift # past argument
      shift # past value
      ;;
    --min_mem)
      min_mem="$2"
      shift # past argument
      shift # past value
      ;;
    --max_mem)
      max_mem="$2"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

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
if [ ! -z "${painting-}" ];
then
	echo "painting   = $painting"
fi
if [ ! -z "${seed-}" ];
then
  echo "seed   = $seed"
fi
echo "P      = ${p}"
echo "q      = ${q}"
echo "min_mem = ${min_mem}"
echo "max_mem = ${max_mem}"
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

######################################################################################################

#### DO NOT CHANGE THE FOLLOWING PARAMETERS
batch_windows=5

######################################################################################################

#### Divide data into chunks of 500k SNPs (the chunks are overlapping to avoid edge effects)
#-sync y causes bsub to wait for the job to complete before exiting.

PATH_TO_RELATE2=${PATH_TO_RELATE}
if [[ "${PATH_TO_RELATE}" != /* ]]; 
then
  PATH_TO_RELATE="../${PATH_TO_RELATE}"
fi
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

    bsub -K -J make_chunks_${output} \
         -env all,PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},dist=${dist},output=${output} \
         -cwd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -M ${min_mem} \
         "${PATH_TO_RELATE}/scripts/RelateLSF/MakeChunks.sh"

  else

    if [[ "$annot" != /* ]]; 
    then
      annot="../${annot}"
    fi

    bsub -K -J make_chunks_${output} \
         -env all,PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},dist=${dist},annot=${annot},output=${output} \
         -cwd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -M ${min_mem} \
         "${PATH_TO_RELATE}/scripts/RelateLSF/MakeChunks.sh"

  fi
else
  if [ -z ${annot-} ]
  then

    bsub -K -J make_chunks_${output} \
         -env all,PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},output=${output} \
         -cwd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -M ${min_mem} \
         "${PATH_TO_RELATE}/scripts/RelateLSF/MakeChunks.sh"

  else

    if [[ "$annot" != /* ]]; 
    then
      annot="../${annot}"
    fi

    bsub -K -J make_chunks_${output} \
         -env all,PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},annot=${annot},output=${output} \
         -cwd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -M ${min_mem} \
         "${PATH_TO_RELATE}/scripts/RelateLSF/MakeChunks.sh"

  fi
fi


###########GET NUMBER OF CHUNKS
num_chunks=($(${PATH_TO_RELATE2}/scripts/RelateLSF/read_bin.py "${output}/${output}/parameters.bin"))
num_chunks=${num_chunks[2]}
prev_chunk=-1

echo "********************************"
echo "Number of chunks: "$num_chunks

######################################################################################################
#### Build trees for each chunk
painting=$(echo ${painting} | awk '{ gsub(",", ";") ; print $0 }'CC)

## paint all sequences against each other
for chunk in `seq 0 $(($num_chunks - 1))`;
do

  #GET NUMBER OF WINDOWS
  num_windows=($(${PATH_TO_RELATE2}/scripts/RelateLSF/read_bin.py "${output}/${output}/parameters_c"${chunk}".bin"))
  num_windows=${num_windows[2]}
  num_batched_windows=$(($num_windows/$batch_windows + 1))

  ## paint all sequences against each other
  #make sure that only 5 paintings exist at a time. Chunk 5 is painted only after chunk 0 is done etc.
 
  if [ ${chunk} -ge 1 ]
	then
		bsub -w find_equivalent_branches_${output}_$(($chunk - 1)) \
				 -J paint_${output}_${chunk} \
				 -cwd ${PWD}/${output} \
				 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=${chunk},painting=${painting},output=${output} \
				 -e log/paint_c${chunk}.log \
				 -o log/paint_c${chunk}.log \
				 -P $p \
				 -q $q \
				 -M ${min_mem} \
       ${PATH_TO_RELATE}/scripts/RelateLSF/Paint.sh
  else
		bsub -J paint_${output}_${chunk} \
				 -cwd ${PWD}/${output} \
				 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=${chunk},output=${output} \
				 -e log/paint_c${chunk}.log \
				 -o log/paint_c${chunk}.log \
				 -P $p \
				 -q $q \
				 -M ${min_mem} \
       ${PATH_TO_RELATE}/scripts/RelateLSF/Paint.sh
	fi
	
	if [ -z ${seed-} ]
	then
		## build tree topologies
		bsub -w paint_${output}_${chunk} \
				 -J build_topology_${output}_${chunk}[1-${num_batched_windows}] \
				 -cwd ${PWD}/${output} \
				 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output},batch_windows=$batch_windows,painting=${painting} \
				 -e build_${output}.log \
				 -o build_${output}.log \
				 -P $p \
				 -q $q \
				 -M ${max_mem} \
				 ${PATH_TO_RELATE}/scripts/RelateLSF/BuildTopology.sh 
  else
		## build tree topologies
		bsub -w paint_${output}_${chunk} \
			-J build_topology_${output}_${chunk}[1-${num_batched_windows}] \
			-cwd ${PWD}/${output} \
			-env all,PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output},batch_windows=$batch_windows,seed=${seed},painting=${painting} \
			-e build_${output}.log \
			-o build_${output}.log \
			-P $p \
			-q $q \
			-M ${max_mem} \
			${PATH_TO_RELATE}/scripts/RelateLSF/BuildTopology.sh 
	fi

  ## find equivalent branches in adjacent trees 
  bsub -w build_topology_${output}_${chunk} \
       -J find_equivalent_branches_${output}_${chunk} \
       -cwd ${PWD}/${output}  \
       -env all,PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output} \
       -e log/find_equivalent_branches_c${chunk}.log \
       -o log/find_equivalent_branches_c${chunk}.log \
       -P $p \
       -q $q \
       -M ${min_mem} \
       ${PATH_TO_RELATE}/scripts/RelateLSF/FindEquivalentBranches.sh

  ## infer branch lengths
  if [ -z ${sample_ages-} ]
	then
		if [ -z ${coal-} ]
		then
			if [ -z ${seed-} ]
			then
				bsub -w find_equivalent_branches_${output}_${chunk} \
						 -J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
						 -cwd ${PWD}/${output} \
						 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -M ${min_mem} \
						 ${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			else
				bsub -w find_equivalent_branches_${output}_${chunk} \
						 -J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
						 -cwd ${PWD}/${output} \
						 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output},seed=${seed} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -M ${min_mem} \
						 ${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			fi

			## combine args into one file
			bsub -w infer_branch_lengths_${output}_${chunk} \
					 -J combine_args_${output} \
					 -cwd ${PWD}/${output} \
					 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,chunk_index=${chunk},output=${output} \
					 -e log/combine_args_c${chunk}.log \
					 -o log/combine_args_c${chunk}.log \
					 -P $p \
					 -q $q \
					 -M ${min_mem} \
					 ${PATH_TO_RELATE}/scripts/RelateLSF/CombineArgs.sh

		else

			if [[ "$coal" != /* ]]; 
			then
				coal="../${coal}"
			fi

			if [ -z ${seed-} ]
			then
				bsub -w find_equivalent_branches_${output}_${chunk} \
						 -J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
						 -cwd ${PWD}/${output} \
						 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -M ${min_mem} \
						 ${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			else
				bsub -w find_equivalent_branches_${output}_${chunk} \
						 -J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
						 -cwd ${PWD}/${output} \
						 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output},seed=${seed} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -M ${min_mem} \
						 ${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			fi

			## combine args into one file
			bsub -w infer_branch_lengths_${output}_${chunk} \
					 -J combine_args_${output} \
					 -cwd ${PWD}/${output} \
					 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},chunk_index=${chunk},output=${output} \
					 -e log/combine_args_c${chunk}.log \
					 -o log/combine_args_c${chunk}.log \
					 -P $p \
					 -q $q \
					 -M ${min_mem} \
					 ${PATH_TO_RELATE}/scripts/RelateLSF/CombineArgs.sh

		fi

  else

		if [ -z ${coal-} ]
		then
			if [ -z ${seed-} ]
			then
				bsub -w find_equivalent_branches_${output}_${chunk} \
					-J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
					-cwd ${PWD}/${output} \
					-env all,PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-M ${min_mem} \
					${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			else
				bsub -w find_equivalent_branches_${output}_${chunk} \
					-J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
					-cwd ${PWD}/${output} \
					-env all,PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output},seed=${seed} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-M ${min_mem} \
					${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			fi

			## combine args into one file
			bsub -w infer_branch_lengths_${output}_${chunk} \
				-J combine_args_${output} \
				-cwd ${PWD}/${output} \
				-env all,PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,chunk_index=${chunk},output=${output} \
				-e log/combine_args_c${chunk}.log \
				-o log/combine_args_c${chunk}.log \
				-P $p \
				-q $q \
				-M ${min_mem} \
				${PATH_TO_RELATE}/scripts/RelateLSF/CombineArgs.sh

		else

			if [[ "$coal" != /* ]]; 
			then
				coal="../${coal}"
			fi

			if [ -z ${seed-} ]
			then
				bsub -w find_equivalent_branches_${output}_${chunk} \
					-J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
					-cwd ${PWD}/${output} \
					-env all,PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-M ${min_mem} \
					${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			else
				bsub -w find_equivalent_branches_${output}_${chunk} \
					-J infer_branch_lengths_${output}_${chunk}[1-${num_batched_windows}] \
					-cwd ${PWD}/${output} \
					-env all,PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output},seed=${seed} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-M ${min_mem} \
					${PATH_TO_RELATE}/scripts/RelateLSF/InferBranchLengths.sh
			fi

			## combine args into one file
			bsub -w infer_branch_lengths_${output}_${chunk} \
				-J combine_args_${output} \
				-cwd ${PWD}/${output} \
				-env all,PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},chunk_index=${chunk},output=${output} \
				-e log/combine_args_c${chunk}.log \
				-o log/combine_args_c${chunk}.log \
				-P $p \
				-q $q \
				-M ${min_mem} \
				${PATH_TO_RELATE}/scripts/RelateLSF/CombineArgs.sh

		fi
	fi

  prev_chunk=$chunk

done

############################# Finalize ###########################
#-sync y causes bsub to wait for the job to complete before exiting. 
#finalize results
if [ -z ${sample_ages-} ]
then
	bsub -K -w combine_args_${output} \
			 -cwd ${PWD}/${output} \
			 -J finalize_${output} \
			 -env all,PATH_TO_RELATE=${PATH_TO_RELATE},output=${output} \
			 -e log/combine_args.log \
			 -o log/combine_args.log \
			 -P $p \
			 -q $q \
			 -M ${min_mem} \
			 ${PATH_TO_RELATE}/scripts/RelateLSF/Finalize.sh
else
	bsub -K -w combine_args_${output} \
		-cwd ${PWD}/${output} \
		-J finalize_${output} \
		-env all,PATH_TO_RELATE=${PATH_TO_RELATE}, sample_ages=${sample_ages},output=${output} \
		-e log/combine_args.log \
		-o log/combine_args.log \
		-P $p \
		-q $q \
		-M ${min_mem} \
		${PATH_TO_RELATE}/scripts/RelateLSF/Finalize.sh
fi

gzip ${output}/${output}.anc
gzip ${output}/${output}.mut

##################################################################
#clean up directory
for chunk in `seq 0 $(($num_chunks - 1))`;
do
  mv ${output}/*.log ${output}/log/
done
cd ${output}
tar -cf log.tar log/
rm -rf log/
cd ..

