#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/RelateSGE/RelateSGE.sh"
  echo ""
  echo "--haps:   Filename of haps file (Output file format of Shapeit)."
  echo "--sample: Filename of sample file (Output file format of Shapeit)."
  echo "--map:    Genetic map."
  echo "-m,--mu:  Mutation rate."
  echo "-N,--Ne:  Effective population size."
  echo "-o,--output: Filename of output without file extension."
  echo "-P:       Assign job to the project."
  echo "-q:       Queue."
  echo "-pe:      Optional. Processor count per node. Default is 'shmem 1'."
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
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/RelateSGE/RelateSGE.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

pe="shmem 1"
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
    --pe)
      pe="$2"
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
echo "pe     = ${pe}"
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
#-sync y causes qsub to wait for the job to complete before exiting.

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

    qsub -sync y \
         -N make_chunks_${output} \
         -v PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},dist=${dist},output=${output} \
         -wd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -pe ${pe} \
         "${PATH_TO_RELATE}/scripts/RelateSGE/MakeChunks.sh"

  else

    if [[ "$annot" != /* ]]; 
    then
      annot="../${annot}"
    fi

    qsub -sync y \
         -N make_chunks_${output} \
         -v PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},dist=${dist},annot=${annot},output=${output} \
         -wd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -pe ${pe} \
         "${PATH_TO_RELATE}/scripts/RelateSGE/MakeChunks.sh"

  fi
else
  if [ -z ${annot-} ]
  then

    qsub -sync y \
         -N make_chunks_${output} \
         -v PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},output=${output} \
         -wd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -pe ${pe} \
         "${PATH_TO_RELATE}/scripts/RelateSGE/MakeChunks.sh"

  else

    if [[ "$annot" != /* ]]; 
    then
      annot="../${annot}"
    fi

    qsub -sync y \
         -N make_chunks_${output} \
         -v PATH_TO_RELATE=${PATH_TO_RELATE},haps=${haps},sample=${sample},map=${map},memory=${memory},annot=${annot},output=${output} \
         -wd ${PWD}/${output} \
         -e log/make_chunks.log \
         -o log/make_chunks.log \
         -P $p \
         -q $q \
         -pe ${pe} \
         "${PATH_TO_RELATE}/scripts/RelateSGE/MakeChunks.sh"

  fi
fi


###########GET NUMBER OF CHUNKS
num_chunks=($(${PATH_TO_RELATE}/scripts/RelateSGE/read_bin.py "${output}/${output}/parameters.bin"))
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
  num_windows=($(${PATH_TO_RELATE}/scripts/RelateSGE/read_bin.py "${output}/${output}/parameters_c"${chunk}".bin"))
  num_windows=${num_windows[2]}
  num_batched_windows=$(($num_windows/$batch_windows + 1))

  ## paint all sequences against each other
  #make sure that only 5 paintings exist at a time. Chunk 5 is painted only after chunk 0 is done etc.
  qsub -hold_jid find_equivalent_branches_${output}_$(($chunk - 5)) \
       -N paint_${output}_${chunk} \
       -wd ${PWD}/${output} \
       -v PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=${chunk},painting=${painting},output=${output} \
       -e log/paint_c${chunk}.log \
       -o log/paint_c${chunk}.log \
       -P $p \
       -q $q \
       -pe ${pe} \
       ${PATH_TO_RELATE}/scripts/RelateSGE/Paint.sh
				 
	if [ -z ${seed-} ]			 
	then
		## build tree topologies
		qsub -hold_jid paint_${output}_${chunk} \
				 -N build_topology_${output}_${chunk} \
				 -wd ${PWD}/${output} \
				 -t 1-$num_batched_windows \
				 -v PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output},batch_windows=$batch_windows,painting=${painting} \
				 -e build_${output}.log \
				 -o build_${output}.log \
				 -P $p \
				 -q $q \
				 -pe ${pe} \
				 ${PATH_TO_RELATE}/scripts/RelateSGE/BuildTopology.sh
  else
		## build tree topologies
		qsub -hold_jid paint_${output}_${chunk} \
			-N build_topology_${output}_${chunk} \
			-wd ${PWD}/${output} \
			-t 1-$num_batched_windows \
			-v PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output},batch_windows=$batch_windows,seed=${seed},painting=${painting} \
			-e build_${output}.log \
			-o build_${output}.log \
			-P $p \
			-q $q \
			-pe ${pe} \
			${PATH_TO_RELATE}/scripts/RelateSGE/BuildTopology.sh
	fi

  ## find equivalent branches in adjacent trees 
  qsub -hold_jid build_topology_${output}_${chunk} \
       -N find_equivalent_branches_${output}_${chunk} \
       -wd ${PWD}/${output}  \
       -v PATH_TO_RELATE=${PATH_TO_RELATE},chunk_index=$chunk,output=${output} \
       -e log/find_equivalent_branches_c${chunk}.log \
       -o log/find_equivalent_branches_c${chunk}.log \
       -P $p \
       -q $q \
       -pe ${pe} \
       ${PATH_TO_RELATE}/scripts/RelateSGE/FindEquivalentBranches.sh

  ## infer branch lengths
  if [ -z ${sample_ages-} ]
	then
		if [ -z ${coal-} ]
		then
			if [ -z ${seed-} ]
			then
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
						 -N infer_branch_lengths_${output}_${chunk} \
						 -wd ${PWD}/${output} \
						 -t 1-$num_batched_windows \
						 -v PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -pe ${pe} \
						 ${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			else
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
						 -N infer_branch_lengths_${output}_${chunk} \
						 -wd ${PWD}/${output} \
						 -t 1-$num_batched_windows \
						 -v PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output},seed=${seed} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -pe ${pe} \
						 ${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			fi

			## combine args into one file
			qsub -hold_jid infer_branch_lengths_${output}_${chunk} \
					 -N combine_args_${output} \
					 -wd ${PWD}/${output} \
					 -v PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,chunk_index=${chunk},output=${output} \
					 -e log/combine_args_c${chunk}.log \
					 -o log/combine_args_c${chunk}.log \
					 -P $p \
					 -q $q \
					 -pe ${pe} \
					 ${PATH_TO_RELATE}/scripts/RelateSGE/CombineArgs.sh

		else

			if [[ "$coal" != /* ]]; 
			then
				coal="../${coal}"
			fi

			if [ -z ${seed-} ]
			then
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
						 -N infer_branch_lengths_${output}_${chunk} \
						 -wd ${PWD}/${output} \
						 -t 1-$num_batched_windows \
						 -v PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -pe ${pe} \
						 ${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			else
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
						 -N infer_branch_lengths_${output}_${chunk} \
						 -wd ${PWD}/${output} \
						 -t 1-$num_batched_windows \
						 -v PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,output=${output},seed=${seed} \
						 -e infer_branch_length_c${chunk}.log \
						 -o infer_branch_length_c${chunk}.log \
						 -P $p \
						 -q $q \
						 -pe ${pe} \
						 ${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			fi

			## combine args into one file
			qsub -hold_jid infer_branch_lengths_${output}_${chunk} \
					 -N combine_args_${output} \
					 -wd ${PWD}/${output} \
					 -v PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},chunk_index=${chunk},output=${output} \
					 -e log/combine_args_c${chunk}.log \
					 -o log/combine_args_c${chunk}.log \
					 -P $p \
					 -q $q \
					 -pe ${pe} \
					 ${PATH_TO_RELATE}/scripts/RelateSGE/CombineArgs.sh

		fi

  else

		if [ -z ${coal-} ]
		then
			if [ -z ${seed-} ]
			then
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
					-N infer_branch_lengths_${output}_${chunk} \
					-wd ${PWD}/${output} \
					-t 1-$num_batched_windows \
					-v PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-pe ${pe} \
					${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			else
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
					-N infer_branch_lengths_${output}_${chunk} \
					-wd ${PWD}/${output} \
					-t 1-$num_batched_windows \
					-v PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output},seed=${seed} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-pe ${pe} \
					${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			fi

			## combine args into one file
			qsub -hold_jid infer_branch_lengths_${output}_${chunk} \
				-N combine_args_${output} \
				-wd ${PWD}/${output} \
				-v PATH_TO_RELATE=${PATH_TO_RELATE},Ne=$Ne,chunk_index=${chunk},output=${output} \
				-e log/combine_args_c${chunk}.log \
				-o log/combine_args_c${chunk}.log \
				-P $p \
				-q $q \
				-pe ${pe} \
				${PATH_TO_RELATE}/scripts/RelateSGE/CombineArgs.sh

		else

			if [[ "$coal" != /* ]]; 
			then
				coal="../${coal}"
			fi

			if [ -z ${seed-} ]
			then
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
					-N infer_branch_lengths_${output}_${chunk} \
					-wd ${PWD}/${output} \
					-t 1-$num_batched_windows \
					-v PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-pe ${pe} \
					${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			else
				qsub -hold_jid find_equivalent_branches_${output}_${chunk} \
					-N infer_branch_lengths_${output}_${chunk} \
					-wd ${PWD}/${output} \
					-t 1-$num_batched_windows \
					-v PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},mu=$mu,chunk_index=$chunk,batch_windows=$batch_windows,sample_ages=${sample_ages},output=${output},seed=${seed} \
					-e infer_branch_length_c${chunk}.log \
					-o infer_branch_length_c${chunk}.log \
					-P $p \
					-q $q \
					-pe ${pe} \
					${PATH_TO_RELATE}/scripts/RelateSGE/InferBranchLengths.sh
			fi

			## combine args into one file
			qsub -hold_jid infer_branch_lengths_${output}_${chunk} \
				-N combine_args_${output} \
				-wd ${PWD}/${output} \
				-v PATH_TO_RELATE=${PATH_TO_RELATE},coal=${coal},chunk_index=${chunk},output=${output} \
				-e log/combine_args_c${chunk}.log \
				-o log/combine_args_c${chunk}.log \
				-P $p \
				-q $q \
				-pe ${pe} \
				${PATH_TO_RELATE}/scripts/RelateSGE/CombineArgs.sh

		fi
	fi

  prev_chunk=$chunk

done

############################# Finalize ###########################
#-sync y causes qsub to wait for the job to complete before exiting. 
#finalize results
if [ -z ${sample_ages-} ]
then
	qsub -sync y \
			 -hold_jid combine_args_${output} \
			 -wd ${PWD}/${output} \
			 -N finalize_${output} \
			 -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output} \
			 -e log/combine_args.log \
			 -o log/combine_args.log \
			 -P $p \
			 -q $q \
			 -pe ${pe} \
			 ${PATH_TO_RELATE}/scripts/RelateSGE/Finalize.sh
else
	qsub -sync y \
		-hold_jid combine_args_${output} \
		-wd ${PWD}/${output} \
		-N finalize_${output} \
		-v PATH_TO_RELATE=${PATH_TO_RELATE}, sample_ages=${sample_ages},output=${output} \
		-e log/combine_args.log \
		-o log/combine_args.log \
		-P $p \
		-q $q \
		-pe ${pe} \
		${PATH_TO_RELATE}/scripts/RelateSGE/Finalize.sh
fi

##################################################################
#clean up directory
mv ${output}/*.log ${output}/log/
cd ${output}
tar -cf log.tar log/
rm -rf log/

gzip ${output}.anc
gzip ${output}.mut

cd ..

